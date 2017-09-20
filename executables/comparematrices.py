import numpy
import sys
import skbio.stats.distance._mantel as mantel

def main():
    if len(sys.argv) != 3:
        print("Usage: python CompareMatrices <matrix1.txt> <matrix2.txt>")
        return
    m1 = numpy.loadtxt(sys.argv[1], dtype=float)
    m2 = numpy.loadtxt(sys.argv[2], dtype=float)
    spr = mantel.mantel(m1, m2, method='spearman')
    prs = mantel.mantel(m1, m2, method='pearson')
    print ("Spearman: ", spr)
    print ("Pearson: ", prs)
    
if __name__ == '__main__':
    main()

