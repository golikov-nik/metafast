package algo;

import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;

import org.apache.log4j.Logger;

import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.executors.NonBlockingQueueExecutor;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.Long2ShortHashMapInterface;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.Tool;
import structures.ConnectedComponent;
import structures.Sequence;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

import static io.IOUtils.withP;

public class ComponentsBuilder {

	public static List<ConnectedComponent> splitStrategy(int alg, BigLong2ShortHashMap hm, Deque<Sequence> sequences,
			int k, int b1, int b2, String statFP, Logger logger, int availableProcessors) throws FileNotFoundException {

		ComponentsBuilder builder = new ComponentsBuilder(k, b1, b2, availableProcessors, statFP, logger);
		switch (alg) {
		case 1:
			builder.runRandomGroupStrategy(hm);
			break;
		case 2:
			builder.runRandomBFSGroupStrategy(hm);
			break;
		case 3:
			builder.runRandomStrategy(hm);
			break;
		case 4:
			builder.runRandomBFSStrategy(hm);
			break;
		case 5:
			builder.runSeptemberAlg(hm, sequences);
			break;
		default:
			builder.run(hm);
			break;
		}
		return builder.ans;
	}

	final private List<ConnectedComponent> ans;
	final private NonBlockingQueueExecutor executor;
	final int k;
	final int b1, b2;
	final String statFP;
	final private Logger logger;

	private ComponentsBuilder(int k, int b1, int b2, int availableProcessors, String statFP, Logger logger) {
		this.ans = new ArrayList<ConnectedComponent>();
		this.executor = new NonBlockingQueueExecutor(availableProcessors, comparator);
		this.k = k;
		this.b1 = b1;
		this.b2 = b2;
		this.statFP = statFP;
		this.logger = logger;
	}

	private void runRandomBFSStrategy(BigLong2ShortHashMap hm) {
		int hmSize = (int) hm.size(); // currently assuming integer number of
										// kmers in graph
		Tool.info(logger, "Kmer groups in graph: " + String.valueOf(hmSize));
		Long[] kmers = new Long[hmSize];
		Iterator<MutableLongShortEntry> it = hm.entryIterator();
		long kmersLeft = 0;
		for (int i = 0; it.hasNext(); i++) {
			MutableLongShortEntry entry = it.next();
			kmersLeft += entry.getValue();
			kmers[i] = entry.getKey();
		}
		Tool.info(logger, "Kmers in graph: " + String.valueOf(kmersLeft));
		Random r = new Random();

		while (kmersLeft > 0) {
			long kmer = kmers[r.nextInt(kmers.length)];
			short value = hm.get(kmer);
			if (value > 0) {
				ConnectedComponent comp = new ConnectedComponent();
				hm.put(kmer, (short) (value - 1));
				Queue<Long> q = new LinkedList<>();
				q.add(kmer);
				while (!q.isEmpty() && comp.size <= b2) {
					long cur = q.poll();
					comp.add(cur);
					kmersLeft--;
					for (long neighbour : neighboursInGraph(hm, cur, k)) {
						short neighbourValue = hm.get(neighbour);
						if (neighbourValue > 0) {
							hm.put(neighbour, (short) (neighbourValue - 1));
							q.add(neighbour);
						}
					}
				}
				for (long notUsedKmer : q) {
					hm.put(notUsedKmer, (short) (hm.get(notUsedKmer) + 1));
				}
				if (comp.size < b1 || comp.size > b2) {
					continue;
				} else {
					Tool.info(logger, "Put component with size: " + comp.size);
					ans.add(comp);
				}
			}
		}
	}

	private void runRandomBFSGroupStrategy(BigLong2ShortHashMap hm) {
		int hmSize = (int) hm.size(); // currently assuming integer number of
										// kmers in graph
		Tool.info(logger, "Kmers in graph: " + String.valueOf(hmSize));
		Random r = new Random();
		HashMap<Long, Boolean> used = new HashMap<>();
		hm.entryIterator().forEachRemaining(e -> used.put(e.getKey(), false));
		Long[] kmers = new Long[hmSize];
		used.keySet().toArray(kmers);
		int kmersLeft = hmSize;
		while (kmersLeft > 0) {
			long kmer = kmers[r.nextInt(kmers.length)];
			if (!used.get(kmer)) {
				ConnectedComponent comp = new ConnectedComponent();
				used.put(kmer, true);
				Queue<Long> q = new LinkedList<>();
				q.add(kmer);
				while (!q.isEmpty() && comp.size <= b2) {
					long cur = q.poll();
					comp.add(cur);
					kmersLeft--;
					for (long neighbour : neighboursInGraph(hm, cur, k)) {
						if (!used.get(neighbour)) {
							used.put(neighbour, true);
							q.add(neighbour);
						}
					}
				}
				for (long notUsedKmer : q) {
					used.put(notUsedKmer, false);
				}
				if (comp.size < b1 || comp.size > b2) {
					continue;
				} else {
					Tool.info(logger, "Put component with size: " + comp.size);
					ans.add(comp);
				}
			}
		}
	}

	private void runSeptemberAlg(BigLong2ShortHashMap hm, Deque<Sequence> sequences) {
		Random r = new Random();
		int hmSize = (int) hm.size();
		HashMap<Long, Boolean> used = new HashMap<>();
		hm.entryIterator().forEachRemaining(e -> used.put(e.getKey(), false));
		Long[] kmers = new Long[hmSize];
		used.keySet().toArray(kmers);
		int kmersLeft = hmSize;
		Tool.info(logger, String.format("Number of sequences: %d\n", sequences.size()));
		int ind = 0;
		int errorKmers = 0;
		for (Sequence seq : sequences) {
			Tool.info(logger, String.format(("Sequence %d length: %d\n"), ind, seq.length() - k + 1));
			ConnectedComponent comp = new ConnectedComponent();
			for (ShortKmer kmer : ShortKmer.kmersOf(seq, k)) {
				if (!used.get(kmer.toLong())) {
					comp.add(kmer.toLong());
					used.put(kmer.toLong(), true);
					kmersLeft--;
				} else {
					errorKmers++;
				}
			}
			if (comp.size > b2) {
				ans.add(comp);
				Tool.info(logger, "Put component with size: " + comp.size);
			}
			ind++;
		}
		Tool.info(logger, String.format("Found %d kmers that were already present", errorKmers));
		while (kmersLeft > 0) {
			long kmer = kmers[r.nextInt(kmers.length)];
			if (!used.get(kmer)) {
				ConnectedComponent comp = new ConnectedComponent();
				used.put(kmer, true);
				Queue<Long> q = new LinkedList<>();
				q.add(kmer);
				while (!q.isEmpty()) {
					long cur = q.poll();
					comp.add(cur);
					kmersLeft--;
					long[] neighbours = neighboursInGraph(hm, cur, k);
					if (comp.size >= b2 && neighbours.length >= 3)
						continue;
					for (long neighbour : neighbours) {
						if (!used.get(neighbour)) {
							used.put(neighbour, true);
							q.add(neighbour);
						}
					}
				}
				for (long notUsedKmer : q) {
					used.put(notUsedKmer, false);
				}
				if (comp.size < b1) {
					continue;
				} else {
					Tool.info(logger, "Put component with size: " + comp.size);
					ans.add(comp);
				}
			}
		}
	}

	public static long[] neighboursInGraph(BigLong2ShortHashMap hm, long kmer, int k) {
		return Arrays.stream(KmerOperations.possibleNeighbours(kmer, k)).filter(e -> hm.get(e) > 0).toArray();
	}

	private void runRandomGroupStrategy(BigLong2ShortHashMap hm) {
		int hmSize = (int) hm.size(); // currently assuming integer number of
										// kmers in graph
		int kmersLeft = hmSize;
		Tool.info(logger, "Kmer groups in graph: " + String.valueOf(hmSize));
		Random r = new Random();
		Long[] kmers = new Long[hmSize];
		Iterator<MutableLongShortEntry> it = hm.entryIterator();
		for (int i = 0; it.hasNext(); i++) {
			MutableLongShortEntry next = it.next();
			kmers[i] = next.getKey();
		}
		ConnectedComponent currentComponent = new ConnectedComponent();
		int currentComponentSize = r.nextInt(Math.min(b2, kmersLeft) - b1 + 1) + b1; // size
		while (kmersLeft > 0) {
			long token = kmers[r.nextInt(hmSize)];
			short value = hm.get(token);
			if (value > 0) {
				if (currentComponent.size == currentComponentSize) {
					ans.add(currentComponent);
					Tool.info(logger, "Put component with size = " + String.valueOf(currentComponent.size));
					if (kmersLeft >= b1) {
						currentComponent = new ConnectedComponent();
						currentComponentSize = r.nextInt(Math.min(b2, kmersLeft) - b1 + 1) + b1;
					} else {
						break;
					}
				}
				currentComponent.add(token);
				hm.put(token, (short) -value);
				kmersLeft--;
			}
		}
		if (currentComponent.size >= b1 && currentComponent.size < currentComponentSize) {
			ans.add(currentComponent);
			Tool.info(logger, "Put component with size = " + String.valueOf(currentComponent.size));
		}
	}

	private void runRandomStrategy(BigLong2ShortHashMap hm) {
		int hmSize = (int) hm.size(); // currently assuming integer number of
										// kmers in graph
		Tool.info(logger, "Kmer groups in graph: " + String.valueOf(hmSize));
		Random r = new Random();
		Long[] kmers = new Long[hmSize];
		Iterator<MutableLongShortEntry> it = hm.entryIterator();
		long sumValue = 0;
		for (int i = 0; it.hasNext(); i++) {
			MutableLongShortEntry next = it.next();
			sumValue += (long) next.getValue();
			kmers[i] = next.getKey();
		}
		Tool.info(logger, "Kmers in graph: " + String.valueOf(sumValue));
		ConnectedComponent currentComponent = new ConnectedComponent();
		assert sumValue >= b1;
		long currentComponentSize = (long) (Math.random() * (Math.min(b2, sumValue) - b1 + 1) + b1); // size
																										// in
																										// range
																										// [b1,
																										// min(b2,
																										// sumValue)]
		while (sumValue > 0) {
			long token = kmers[r.nextInt(hmSize)];
			short value = hm.get(token);
			if (value > 0) {
				if (currentComponent.size == currentComponentSize) {
					ans.add(currentComponent);
					Tool.info(logger, "Put component with size = " + String.valueOf(currentComponent.size));
					if (sumValue >= b1) {
						currentComponent = new ConnectedComponent();
						currentComponentSize = (long) (Math.random() * (Math.min(b2, sumValue) - b1 + 1) + b1);
					} else {
						break;
					}
				}
				currentComponent.add(token);
				hm.put(token, (short) (value - 1));
				sumValue--;
			}
		}
		if (currentComponent.size >= b1 && currentComponent.size < currentComponentSize) {
			ans.add(currentComponent);
			Tool.info(logger, "Put component with size = " + String.valueOf(currentComponent.size));
		}
	}

	private void run(BigLong2ShortHashMap hm) throws FileNotFoundException {
		Tool.info(logger, "First iteration...");
		Timer t = new Timer();
		long hmSize = hm.size();
		int curFreqThreshold = 1; // current component is formed of k-mers with
									// frequency >= 1
		List<ConnectedComponent> newComps = findAllComponents(hm, k, b2, curFreqThreshold);

		int small = 0, ok = 0, big = 0;
		long smallK = 0, okK = 0;

		List<ConnectedComponent> toProcess = new ArrayList<ConnectedComponent>();
		for (ConnectedComponent comp : newComps) {
			if (comp.size < b1) {
				small++;
				smallK += comp.size;
			} else if (comp.size <= b2) {
				ok++;
				okK += comp.size;
				ans.add(comp);
			} else {
				big++;
				toProcess.add(comp);
			}
		}

		int ansFirst = ans.size();
		Tool.info(logger, "Found " + NumUtils.groupDigits(ok) + " good components, " + "and "
				+ NumUtils.groupDigits(big) + " big ones");
		Tool.info(logger, "First iteration was finished in " + t);

		Tool.debug(logger, "Total components found = " + NumUtils.groupDigits(newComps.size()) + ", " + "kmers = "
				+ NumUtils.groupDigits(hmSize));
		Tool.debug(logger, "Components count: small = " + withP(small, newComps.size()) + ", " + "ok = "
				+ withP(ok, newComps.size()) + ", " + "big = " + withP(big, newComps.size()));
		Tool.debug(logger, "Components kmers: small = " + withP(smallK, hmSize) + ", " + "ok = " + withP(okK, hmSize)
				+ ", " + "big = " + withP(hmSize - smallK - okK, hmSize));
		Tool.debug(logger, "FreqThreshold = " + curFreqThreshold + ", " + "components added = " + ok
				+ ", total components added = " + ans.size());

		Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", "
				+ "after it = " + Misc.usedMemoryAsString());

		hm = null; // for cleaning
		newComps = null;

		Tool.debug(logger, "Memory used after cleaning = " + Misc.usedMemoryAsString() + ", final time = " + t);

		if (big != 0) {
			Tool.info(logger, "Following iterations...");
			t.start();

			for (ConnectedComponent comp : toProcess) {
				executor.addTask(new Task(comp));
			}
			toProcess = null;

			ConnectedComponent biggest = ((Task) executor.tasks.peek()).component;
			Tool.debug(logger, "Biggest component has " + withP(biggest.size, hmSize, "kmers", "of initial hm size"));
			Tool.debug(logger,
					"Saved to new hm from it = " + withP(biggest.nextHM.size(), biggest.size, "kmers", "of its size"));
			biggest = null;

			executor.startWorkers();

			executor.waitForTasksToFinish();
			executor.shutdownAndAwaitTermination();

			Tool.info(logger, "Found " + NumUtils.groupDigits(ans.size() - ansFirst) + " good components "
					+ "extracted from big components");
			Tool.info(logger, "All the following iterations were finished in " + t);
			Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", "
					+ "after it = " + Misc.usedMemoryAsString());
		}

		// post processing...
		Tool.debug(logger, "ans.size = " + ans.size());

		Collections.sort(ans);

		PrintWriter statPW = new PrintWriter(statFP);
		statPW.println("# component.no\tcomponent.size\tcomponent.weight\tusedFreqThreshold");
		for (int i = 0; i < ans.size(); i++) {
			ConnectedComponent comp = ans.get(i);
			statPW.println((i + 1) + "\t" + comp.size + "\t" + comp.weight + "\t" + comp.usedFreqThreshold);
		}
		statPW.close();
	}

	class Task implements Runnable {
		ConnectedComponent component; // task is to split this component to good
										// ones

		Task(ConnectedComponent component) {
			this.component = component;
		}

		@Override
		public void run() {
			int curFreqThreshold = component.usedFreqThreshold + 1;

			List<ConnectedComponent> newComps = findAllComponents(component.nextHM, k, b2, curFreqThreshold);

			for (ConnectedComponent comp : newComps) {
				if (comp.size < b1) {
					// skipping
				} else if (comp.size <= b2) {
					synchronized (ans) {
						ans.add(comp);
					}
				} else {
					executor.addTask(new Task(comp));
				}
			}
		}
	}

	Comparator<Runnable> comparator = new Comparator<Runnable>() {
		@Override
		public int compare(Runnable o1, Runnable o2) {
			if (!(o1 instanceof Task) || !(o2 instanceof Task)) {
				return 0;
			}
			Task t1 = (Task) o1;
			Task t2 = (Task) o2;
			return -Long.compare(t1.component.size, t2.component.size);
		}
	};

	/**
	 * Assuming running in one thread for current hm!
	 */
	private static List<ConnectedComponent> findAllComponents(Long2ShortHashMapInterface hm, int k, int b2,
			int curFreqThreshold) {
		List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();
		LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size() / 2));

		Iterator<MutableLongShortEntry> iterator = hm.entryIterator();
		while (iterator.hasNext()) {
			MutableLongShortEntry startKmer = iterator.next();
			if (startKmer.getValue() > 0) { // i.e. if not precessed
				ConnectedComponent comp = bfs(hm, startKmer.getKey(), queue, k, b2, curFreqThreshold);
				ans.add(comp);
			}
		}

		return ans;
	}

	/**
	 * Breadth-first search to make the traversal of the component. If the
	 * component is small (less than b2 vertices), all its kmers is saved to
	 * ConnectedComponent.kmers, else a subset of hm is stored to
	 * ConnectedComponent.nextHM structure.
	 */
	private static ConnectedComponent bfs(Long2ShortHashMapInterface hm, long startKmer, LongArrayFIFOQueue queue,
			int k, int b2, int curFreqThreshold) {

		ConnectedComponent comp = new ConnectedComponent();
		comp.usedFreqThreshold = curFreqThreshold;

		queue.clear();

		queue.enqueue(startKmer);
		short value = hm.get(startKmer);
		assert value > 0;
		hm.put(startKmer, (short) -value); // removing
		comp.add(startKmer, value);
		boolean alreadyBigComp = false;

		while (queue.size() > 0) {
			long kmer = queue.dequeue();

			for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
				value = hm.get(neighbour);
				if (value > 0) { // i.e. if not precessed
					queue.enqueue(neighbour);
					hm.put(neighbour, (short) -value);

					if (!alreadyBigComp) {
						comp.add(neighbour, value);
						if (comp.size > b2) {
							alreadyBigComp = true;
							comp.nextHM = new BigLong2ShortHashMap(4, 13);
							for (long kk : comp.kmers) {
								value = (short) -hm.get(kk);
								assert value > 0;
								if (value >= curFreqThreshold + 1) {
									comp.nextHM.put(kk, value);
								}
							}
							comp.kmers = null;
						}
					} else {
						if (value >= curFreqThreshold + 1) { // for next HM
							comp.nextHM.put(neighbour, value);
						}
						comp.size++;
					}
				}
			}
		}

		return comp;
	}

}
