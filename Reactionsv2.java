import java.io.*;
import java.util.*;
/***
Version 2 of Reactions, trying to clean up code so debugging is easier
*/
public class Reactionsv2 {
    public Molecule reactants;
    public Molecule products;

    public static void main(String[] args) throws FileNotFoundException {
	Reactionsv2 reactions = new Reactionsv2();
	reactions.createReaction(args);
	Molecule reactant = reactions.reactants;
	for (int i = 0; i < reactant.getAtoms().size(); i++) {
	    System.out.print(reactant.getAtoms().get(i) + " ");
	}
	System.out.println();
	reactant.generateFreq();
	int [] freq = reactant.getFreqs();
	for (int i=0; i<freq.length; i++) {
	    System.out.print(freq[i] + " ");
	}
	System.out.println();
	System.out.println("molecules = " + reactant.findNumOfMlcs());
	HashMap<String, Integer> bondTypes = reactant.listChangedBonds();
	Set<String> bondPics = bondTypes.keySet();
	for (String s : bondPics) {
	    System.out.println(s + " " + bondTypes.get(s));
	}
	reactant.buildMatchTrees();
	int atomNumCount = 1;
	for (ArrayList<String> atomTree : reactant.getMatchTrees()) {
	    System.out.print(atomNumCount + ": ");
	    for (String s : atomTree) {
		System.out.print(s + ", ");
	    }
	    System.out.println();
	    atomNumCount++;
	}
	Molecule product = reactions.products;
	product.generateFreq();
	int [] pfreq = product.getFreqs();
	for (int i=0; i<pfreq.length; i++) {
	    System.out.print(pfreq[i] + " ");
	}
	System.out.println();
	System.out.println("molecules = " + product.findNumOfMlcs());
	HashMap<String, Integer> pbondTypes = product.listChangedBonds();
	Set<String> pbondPics = pbondTypes.keySet();
	for (String s : pbondPics) {
	    System.out.println(s + " " + pbondTypes.get(s));
	}
	product.buildMatchTrees();
	int patomNumCount = 1;
	for (ArrayList<String> patomTree : product.getMatchTrees()) {
	    System.out.print(patomNumCount + ": ");
	    for (String s : patomTree) {
		System.out.print(s + ", ");
	    }
	    System.out.println();
	    patomNumCount ++;
	}
	reactions.matchAtoms(reactant, product);
    }
    
    public void createReaction(String[] args) throws FileNotFoundException {
	File reactant;
	Scanner reactantSc;
	Molecule rmolecule;
	File product;
	Scanner productSc;
	Molecule pmolecule;
	if (args.length != 2) {
	    System.out.println("Wrong number  of arguments, enter 2 arguments: reactant product");
	    System.exit(0);
	}
	reactant = new File (args[0]);
	reactantSc = new Scanner(reactant);
	rmolecule = createMolecule(reactantSc);
	reactants = rmolecule;
	product = new File(args[1]);
	productSc = new Scanner(product);
	pmolecule = createMolecule(productSc);
	products = pmolecule;
	System.out.println("Reactant # of atoms = " + reactants.getNumberOfAtoms());
	System.out.println("Reactant # of bonds = " + reactants.getNumberOfBonds());
	System.out.println("Reactant atoms      = " + reactants.getAtoms());
	System.out.println("Reactant bonds      = " + reactants.getBonds());
	System.out.println("Product # of atoms = " + products.getNumberOfAtoms());
	System.out.println("Product # of bonds = " + products.getNumberOfBonds());
	System.out.println("Product atoms      = " + products.getAtoms());
	System.out.println("Product bonds      = " + products.getBonds());
    }
    
    public Molecule createMolecule(Scanner molsc) {
	ArrayList<String> bonds = new ArrayList<String>();
	ArrayList<String> atoms = new ArrayList<String>();
	String info = "";
	String trash = molsc.nextLine();
	Scanner trashsc = new Scanner(trash);
	Boolean molFileStart = false;
	while (molFileStart == false) {
	    if (trashsc.hasNextInt()) {
		info = trash;
		molFileStart = true;
	    } else {
		trash = molsc.nextLine();
		trashsc = new Scanner(trash);
	    }
	}
	Scanner infosc = new Scanner (info);
	int numberOfAtoms = infosc.nextInt();
	int numberOfBonds = infosc.nextInt();
	String temp = molsc.nextLine();
	while (temp.charAt(2) == ' ') {
	    Scanner tempsc = new Scanner(temp);
	    String atom = tempsc.next();
	    atom = tempsc.next();
	    atom = tempsc.next();
	    atom = tempsc.next();
	    atoms.add(atom);
	    temp = molsc.nextLine();
	}
	Scanner bondsc = new Scanner(temp);
	while (bondsc.hasNextInt()) {
	    String bond = "" + bondsc.nextInt() + " " + bondsc.nextInt() + " " + bondsc.nextInt();
	    bonds.add(bond);
	    temp = molsc.nextLine();
	    bondsc = new Scanner(temp);
	}
	Molecule molecule = new Molecule(atoms, numberOfBonds, numberOfAtoms, bonds);
	return molecule;
    }

    public void matchAtoms(Molecule m1, Molecule m2) {
	Map <Integer, Integer> mapping = new HashMap<Integer, Integer>();
	int m1AtomNum = 0;
	int m2AtomNum = 0;
	int rootMatches = 0;
	HashMap<String, Integer> m1BondTypes = m1.listChangedBonds();
	HashMap<String, Integer> m2BondTypes = m2.listChangedBonds();
	Set<String> m1Keys = m1BondTypes.keySet();
	HashMap<String, ArrayList<String>> changedBonds = new HashMap<String, ArrayList<String>>();
	changedBonds.put("+", new ArrayList<String>());
	changedBonds.put("-", new ArrayList<String>());
	boolean firstThrough = false;
	for (String key : m1Keys) {
	    if (!m2BondTypes.containsKey(key)) {
		changedBonds.get("-").add(key);
	    } else {
		for (String key2 : m2BondTypes.keySet()) {
		    if (!m1BondTypes.containsKey(key2) && firstThrough == false) {
			changedBonds.get("+").add(key2);
		    } else if (!m1BondTypes.containsKey(key2) && firstThrough == true) {
			continue;
		    } else {
			int key1Freq = m1BondTypes.get(key);
			int key2Freq = m2BondTypes.get(key);
			if (key1Freq - key2Freq > 0) {
			    changedBonds.get("-").add(key);
			} else if (key1Freq - key2Freq < 0) {
			    changedBonds.get("+").add(key);
			}
		    }
		} // close key2 for loop
	    }
	} // close key for loop
	ArrayList<String> rmvedNodes = new ArrayList<String>();
	ArrayList<String> addedNodes = new ArrayList<String>();
	for (String add : changedBonds.get("+")) {
	    //System.out.println("+ " + add);
	    String parentAtom = "" + add.charAt(0);
	    String bondedAtom = "" + add.charAt(2);
	    String bondtype = "" + add.charAt(1);
	    if (bondtype.equals("-")) {
		addedNodes.add(bondedAtom + "1");
	    } else if (bondtype.equals("=")) {
		addedNodes.add(bondedAtom + "2");
	    } else {
		addedNodes.add(bondedAtom + "3");
	    }
	} // close add for loop
	for (String lose : changedBonds.get("-")) {
	    //System.out.println("- " + lose);
	    String  parentAtom = "" + lose.charAt(0);
	    String bondedAtom = "" + lose.charAt(2);
	    String bondtype = "" + lose.charAt(1);
	    if (bondtype.equals("-")) {
	        rmvedNodes.add(bondedAtom + "1");
	    } else if (bondtype.equals("=")) {
		rmvedNodes.add(bondedAtom + "2");
	    } else {
		rmvedNodes.add(bondedAtom + "3");
	    }
	} // close lose for loop
	System.out.println("Added Nodes: " + addedNodes);
	System.out.println("Lost Nodes: " + rmvedNodes);
	
	ArrayList<Integer> m1atomsMapped = new ArrayList<Integer>();
	ArrayList<Integer> m2atomsMapped = new ArrayList<Integer>();
	HashMap<Integer, Integer> atomMap = new HashMap<Integer, Integer>();
	ArrayList<matchTreeNode> m1trees = m1.getatomTrees();
	ArrayList<matchTreeNode> m2trees = m2.getatomTrees();
	// match same trees first
	for (int trackm1 = 0; trackm1 < m1trees.size(); trackm1++) {
	    for (int trackm2 = 0; trackm2 < m2trees.size(); trackm2++) {
		matchTreeNode m1branch = m1trees.get(trackm1);
		matchTreeNode m2branch = m2trees.get(trackm2);
		if (m2atomsMapped.contains(trackm2)) {
		    continue;
		} else {
		    boolean match = treeMatch(m1branch, m2branch);
		    if (match == true) {
			m1atomsMapped.add(trackm1);
			m2atomsMapped.add(trackm2);
		    }
		}
	    }
	}
    } // close matchAtoms method

    public boolean treeMatch(matchTreeNode m1, matchTreeNode m2) {
	boolean match = false;
	if (m1.equals(m2)) {
	    ArrayList<matchTreeNode> m1level2 = m1.getChildren();
	    ArrayList<matchTreeNode> m2level2 = m2.getChildren();
	    if (m1level2.size() != m2level2.size()) {
		return false;
	    }
	    match = compareChildren(m1level2, m2level2);
	    if (match == false){
		return false;
	    }
	    ArrayList<Integer> m2lvl2matched = new ArrayList<Integer>();
	    for (matchTreeNode m1lvl2 : m1level2) {
		while (match != false) {
		    for (int m2lvl2track = 0; m2lvl2track < m2level2.size(); m2lvl2track++) {
			matchTreeNode m2lvl2 = m2level2.get(m2lvl2track);
			if (m1lvl2.toString().equals(m2lvl2.toString())) {
			    if (compareChildren(m1lvl2.getChildren(), m2lvl2.getChildren()) == true && !m2lvl2matched.contains(m2lvl2)){
				m2lvl2matched.add(m2lvl2track);
				ArrayList<matchTreeNode> m1level3 = m1lvl2.getChildren();
				ArrayList<matchTreeNode> m2level3 = m2lvl2.getChildren();
				match = compareChildren(m1level3, m2level3);
				if (match == false) {
				    return false;
				} else {
				    ArrayList<Integer> m2lvl3matched = new ArrayList<Integer>();
				    
				}
				break;
			    } else if(m2lvl2track+1 != m2level2.size())  {
				continue;
			    } else {
				match = false;
			    }
			}
		    }
		}
		if (match == false) {
		    return false;
		}
	    }
	    
	    
	    
	} else {
	    return false;
	}
	return match;
    } // close treeMatch method

    public boolean compareChildren(ArrayList<matchTreeNode> m1, ArrayList<matchTreeNode> m2) {
	ArrayList<matchTreeNode> unMatchedNodes = new ArrayList<matchTreeNode>();
	ArrayList<Integer> m2visited = new ArrayList<Integer>();
	if (m1.size() != m2.size()) {
	    return false;
	}
	for (int i=0; i<m1.size(); i++) {
	    String temp = m1.get(i).toString();
	    for (int j=0; j<m2.size(); j++) {
		if (temp.equals(m2.get(j).toString()) && !m2visited.contains(j)) {
		    m2visited.add(j);
		    break;
		} else {
		    continue;
		}
	    }
	}
	if (m2visited.size() == m2.size()) {
	    return true;
	} else {
	    return false;
	}
    } // close compareChildren method
}