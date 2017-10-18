import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import structures.ConnectedComponent;

@SuppressWarnings("serial")
public class SimpleGUI extends JFrame implements MouseWheelListener, MouseMotionListener {
	static ArrayList<ArrayList<Kmer>> genomes;
	static ArrayList<Integer> genomeSizes = new ArrayList<>();
	static JScrollPane jsp;
	static int size = 0;
	static int defsize = 0;
	static Color[] clrs;
	static int genomeNumber = 0;
	static Color BGcolor = Color.BLACK;
	static List<ConnectedComponent> components;
	static JPanel drawpanel = new JPanel() {

		@Override
		public void paint(java.awt.Graphics g) {
			BufferedImage image = new BufferedImage(drawpanel.getWidth(), drawpanel.getHeight(),
					BufferedImage.TYPE_INT_RGB);
			Graphics g2 = image.getGraphics();
			int rcsz = this.getWidth() - 60; // rectangle size
			g2.setColor(BGcolor);
			g2.fillRect(0, 0, drawpanel.getWidth(), drawpanel.getHeight());
			g2.setColor(Color.gray);
			for (int i = 0; i < genomeNumber; i++) {
				g2.fillRect(30, 30 * (2 * i + 1), drawpanel.getWidth() - 60, 30);
			}
			int ind = 0;
			for (ArrayList<Kmer> genome : genomes) {
				for (Kmer kmer : genome) {
					g2.setColor(kmer.color);
					double coef = ((double) kmer.position / (double) genomeSizes.get(ind));
					g2.fillRect(30 + (int) (coef * rcsz), 30 * (2 * ind + 1), 1, 30);
				}
				ind++;
			}
			g.drawImage(image, 0, 0, null);
		}
	};

	static class Kmer {
		int position;
		Color color;

		Kmer(int position, Color color) {
			this.position = position;
			this.color = color;
		}

	}

	public SimpleGUI() {
		super("TEST01");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setExtendedState(JFrame.MAXIMIZED_BOTH);

		drawpanel.setLayout(null);
		addMouseWheelListener(this);
		addMouseMotionListener(this);
		drawpanel.setBounds(0, 0, Frame.MAXIMIZED_HORIZ + size, 60 * (genomeNumber + 1));
		if (Frame.MAXIMIZED_VERT > 60 * (genomeNumber + 1)) {
			drawpanel.setBounds(0, 0, Frame.MAXIMIZED_HORIZ + size, Frame.MAXIMIZED_VERT);
		}
		drawpanel.setVisible(true);
		drawpanel.setFocusable(true);
		drawpanel.setPreferredSize(new Dimension(Frame.MAXIMIZED_HORIZ + size, 60 * (genomeNumber + 1)));
		jsp = new JScrollPane(drawpanel);
		jsp.addMouseWheelListener(this);
		jsp.addMouseMotionListener(this);
		jsp.setPreferredSize(new Dimension(Frame.MAXIMIZED_HORIZ + size, 60 * (genomeNumber + 1)));
		if (Frame.MAXIMIZED_VERT > 60 * (genomeNumber + 1)) {
			drawpanel.setPreferredSize(new Dimension(Frame.MAXIMIZED_HORIZ + size, Frame.MAXIMIZED_VERT));
			drawpanel.setBounds(0, 0, Frame.MAXIMIZED_HORIZ + size, Frame.MAXIMIZED_VERT);
			jsp.setPreferredSize(new Dimension(Frame.MAXIMIZED_HORIZ + size, Frame.MAXIMIZED_VERT));
		}

		add(jsp);
	}

	public static void main(String[] args) throws ExecutionFailedException, IOException {
		if (args.length != 3) {
			System.out.println("Usage: Selector <genomes folder> <components.bin> <k>");
			return;
		}
		int k = Integer.parseInt(args[2]);
		File genomesFolder = new File(args[0]);
		assert genomesFolder.isDirectory();
		components = ConnectedComponent.loadComponents(new File(args[1]));
		System.out.printf("Loaded %d components\n\n", components.size());
		Random r = new Random();
		clrs = new Color[components.size()];
		for (int i = 0; i < components.size(); i++) {
			clrs[i] = new Color(r.nextInt(255), r.nextInt(255), r.nextInt(255));
		}
		genomes = new ArrayList<ArrayList<Kmer>>();
		for (File genFile : genomesFolder.listFiles()) {
			if (isValidInputFilename(genFile.getName())) {
				ArrayList<Kmer> genome = new ArrayList<>();
				int ind = 0;
				List<Dna> reads = ReadersUtils.loadDnas(genFile);
				for (Dna dna : reads) {
					for (ShortKmer skmer : ShortKmer.kmersOf(dna, k)) {
						long val = skmer.toLong();
						int comp_ind = 0;
						for (ConnectedComponent comp : components) {
							if (comp.kmers.contains(val)) {
								Kmer kmer = new Kmer(ind, clrs[comp_ind]);
								genome.add(kmer);
							}
							comp_ind++;
						}
						ind++;
					}
				}
				genomeSizes.add(ind + 1);
				genomes.add(genome);
				genomeNumber++;
				System.out.printf("Loaded genome with size=%d\n", ind + 1);
				System.out.printf("Genome has %d kmers found in any of components\n\n", genome.size());
			}
		}
		System.out.println();
		System.out.println("Starting drawing...");
		SimpleGUI app = new SimpleGUI();
		app.setVisible(true);

	}

	public static boolean isValidInputFilename(String filename) {
		return filename.endsWith(".fa") || filename.endsWith(".fasta") || filename.endsWith(".fastq")
				|| filename.endsWith(".fna");
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		if (e.getWheelRotation() > 0) {
			size += 10;
		} else {
			size -= 10;
		}
		if (defsize == 0) {
			defsize = drawpanel.getWidth();
		}
		drawpanel.setPreferredSize(new Dimension(defsize + size, drawpanel.getHeight()));
		drawpanel.setBounds(0, 0, defsize + size, drawpanel.getHeight());
		jsp.setFocusable(false);
		jsp.setPreferredSize(new Dimension(defsize + size, drawpanel.getHeight()));
		drawpanel.repaint();
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseMoved(MouseEvent e) {

	}

}