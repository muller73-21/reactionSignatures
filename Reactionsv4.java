import java.io.*;
import java.util.*;
/***
Version 4 of Reactions, trying to clean up code so debugging is easier
*/
public class Reactionsv4 {
    public Molecule reactants;
    public Molecule products;

    public static void main(String[] args) throws FileNotFoundException {
	Reactionsv4 reactions = new Reactionsv4();
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
	    Scanner bondsc = new Scanner(add);
	    String parentAtom = bondsc.next();
	    String bondtype = bondsc.next();
	    String bondedAtom = bondsc.next();	   
	    if (bondtype.equals("-")) {
		addedNodes.add(bondedAtom + "1");
		addedNodes.add(parentAtom + "1");
	    } else if (bondtype.equals("=")) {
		addedNodes.add(bondedAtom + "2");
	    } else {
		addedNodes.add(bondedAtom + "3");
	    }
	} // close add for loop
	for (String lose : changedBonds.get("-")) {
	    //System.out.println("- " + lose);
	    Scanner bondsc = new Scanner(lose);
	    String  parentAtom = bondsc.next();
	    String bondtype = bondsc.next();
	    String bondedAtom = bondsc.next();	    
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
	
       	HashMap<Integer, Integer> atomMap = new HashMap<Integer, Integer>();
	ArrayList<matchTreeNode> m1parents = m1.getatomTrees();
	ArrayList<matchTreeNode> m2parents = m2.getatomTrees();
	ArrayList<Integer> m2visited = new ArrayList<Integer>();
	ArrayList<Integer> m2eqMatch = new ArrayList<Integer>();
	ArrayList<Integer> m1eqMatch = new ArrayList<Integer>();
	for (int m1lvl1 = 0; m1lvl1 < m1parents.size(); m1lvl1++) {
	    m2loop:
	    for (int m2lvl1 = 0; m2lvl1 < m2parents.size(); m2lvl1++) {
		matchTreeNode m1root = m1parents.get(m1lvl1);
		matchTreeNode m2root = m2parents.get(m2lvl1);
		if (m1root.toString().equals(m2root.toString()) && !m2eqMatch.contains(m2lvl1)) {
		    ArrayList<matchTreeNode> m1children = m1root.getChildren();
		    ArrayList<matchTreeNode> m2children = m2root.getChildren();
		    //System.out.println(m1lvl1 + " " + m2lvl1 + " " + m1children + " " + m2children);
		    ArrayList<matchTreeNode> childrenTest = compareChildren(m1children, m2children);
		    boolean ctest = childrenTest.size() == 0;
		    //System.out.println("children test = " + ctest);
		    if (childrenTest.size() == 0) {
			ArrayList<Integer> m2childrenMatched = new ArrayList<Integer>();
			for (int m1lvl2 = 0; m1lvl2 < m1children.size(); m1lvl2++) {
			    for (int m2lvl2 = 0; m2lvl2 < m2children.size(); m2lvl2++) {
				matchTreeNode m1lvl2curr = m1children.get(m1lvl2);
				matchTreeNode m2lvl2curr = m2children.get(m2lvl2);
				if (m1lvl2curr.toString().equals(m2lvl2curr.toString()) && !m2childrenMatched.contains(m2lvl2)) {
				    ArrayList<matchTreeNode> m1gchildren = m1lvl2curr.getChildren();
				    ArrayList<matchTreeNode> m2gchildren = m2lvl2curr.getChildren();
				    //System.out.println(m1lvl1 + " " + m2lvl1 + " " + m1gchildren + " " + m2gchildren);
				    ArrayList<matchTreeNode> gchildrenTest = compareChildren(m1gchildren, m2gchildren);
				    if (gchildrenTest.size() == 0) {
					if (m1gchildren.size() == 0 && m2gchildren.size() == 0) {
					    m2childrenMatched.add(m2lvl2);
					    //System.out.println(m1lvl1 + " " + m2lvl1 + " " + m2childrenMatched.size() + " " + m2children.size());
					    if (m2childrenMatched.size() == m2children.size()) {
						int i = m1lvl1 + 1;
						int i2 = m2lvl1 + 1;
						mapping.put(i, i2);
						m1eqMatch.add(m1lvl1);
						m2eqMatch.add(m2lvl1);
						break m2loop;
					    }
					}
					ArrayList<Integer> m2gchildrenMatched = new ArrayList<Integer>();
					for (int m1lvl3 = 0; m1lvl3 < m1gchildren.size(); m1lvl3 ++) {
					    for (int m2lvl3 = 0; m2lvl3 < m2gchildren.size(); m2lvl3++) {
						matchTreeNode m1lvl3curr = m1gchildren.get(m1lvl3);
						matchTreeNode m2lvl3curr = m2gchildren.get(m2lvl3);
						if (m1lvl3curr.toString().equals(m2lvl3curr.toString()) && !m2gchildrenMatched.contains(m2lvl3)) {
						    ArrayList<matchTreeNode> m1ggchildren = m1lvl3curr.getChildren();
						    ArrayList<matchTreeNode> m2ggchildren = m2lvl3curr.getChildren();
						    ArrayList<matchTreeNode> ggchildrenTest = compareChildren(m1ggchildren, m2ggchildren);
						    if (ggchildrenTest.size() == 0) {
							m2gchildrenMatched.add(m2lvl3);
							if (m2gchildrenMatched.size() == m2gchildren.size()) {
							    m2childrenMatched.add(m2lvl2);
							    if (m2childrenMatched.size() == m2children.size()) {
								int index = m1lvl1 + 1;
								int index2 = m2lvl1 + 1;
								mapping.put(index, index2);
								m1eqMatch.add(m1lvl1);
								m2eqMatch.add(m2lvl1);
								break m2loop;
							    }							    
							}
						    } 
						} 
					    }
					} 
				    }
				} 
			    }
			}
		    } 
		}
	    }
	} // eq test
	/*for (int m1atom : mapping.keySet()) {
	    System.out.println("reactant atom " + m1atom + " maps to " + mapping.get(m1atom) + " in the product");
	}
	*/
	if (addedNodes.size() == 0 && rmvedNodes.size() == 0) {
	    // match testing
	    for (int m1parIndex = 0; m1parIndex < m1parents.size(); m1parIndex++) {
		matchTreeNode reactantTest = m1parents.get(m1parIndex);
		ArrayList<Integer> m1bonds = reactantTest.getAtomBonds();
		int atomNum = reactantTest.getAtomNumber();
		//System.out.println("Reactant atom " + atomNum + " is bonded to " + m1bonds);
		ArrayList<Integer> possible = new ArrayList<Integer>();
		for (int i=0; i<m1bonds.size(); i++){
		    if (mapping.keySet().contains(m1bonds.get(i)) && !mapping.keySet().contains(atomNum)) {
			int bondedAtom = mapping.get(m1bonds.get(i));
			int atomNum2 = m2parents.get(bondedAtom).getAtomNumber()-1;
			//System.out.println("Product atom " + atomNum2 + " is bonded to " + m2parents.get(atomNum2-1).getAtomBonds() + " is from " + m1bonds.get(i));
			ArrayList<Integer> m2bonds = m2parents.get(atomNum2-1).getAtomBonds();
			for (Integer k : m2bonds) {
			    possible.add(k);
			}
			for (int y = 0; y < possible.size(); y ++) {
			    int j = possible.get(y);
			    if (m2eqMatch.contains(j-1)) {
				//System.out.println(j-1 + " " + m2eqMatch);
				possible.remove(y);
				y --;
			    }
			}
			if (possible.size() == 1) {
			    //System.out.println("match");
			    mapping.put(atomNum, possible.get(0));
			    m1eqMatch.add(atomNum-1);
			    m2eqMatch.add(possible.get(0)-1);
			} if (possible.size() > 1) {
			    for (int q = 0; q < possible.size(); q++) {
				//System.out.println(reactantTest + " " + m2parents.get(possible.get(q)-1) + " " + possible.get(q));
				if (reactantTest.toString().equals(m2parents.get(possible.get(q)-1).toString())) {
				    mapping.put(atomNum, possible.get(q));
				    m1eqMatch.add(atomNum - 1);
				    m2eqMatch.add(possible.get(q)-1);
				}
			    }
			}
		    }
		}
	    }
	    /*
	      System.out.println();
	      for (int m1atom : mapping.keySet()) {
	      System.out.println("reactant atom " + m1atom + " maps to " + mapping.get(m1atom) + " in the product");
	      }
	      System.out.println(m1eqMatch);
	      System.out.println(m2eqMatch);
	    */
	    for (int m1parIndex = 0; m1parIndex < m1parents.size(); m1parIndex++) {
		matchTreeNode reactantTest = m1parents.get(m1parIndex);
		ArrayList<Integer> m1bonds = reactantTest.getAtomBonds();
		int atomNum = reactantTest.getAtomNumber();
		
		//System.out.println("Reactant atom " + atomNum + " is bonded to " + m1bonds);
		ArrayList<Integer> possible = new ArrayList<Integer>();
		for (int i=0; i<m1bonds.size(); i++){
		    
		    if (mapping.keySet().contains(m1bonds.get(i)) && !mapping.keySet().contains(atomNum)) {
			int bondedAtom = mapping.get(m1bonds.get(i));
			int atomNum2 = m2parents.get(bondedAtom).getAtomNumber()-1;
			//System.out.println("Product atom " + atomNum2 + " is bonded to " + m2parents.get(atomNum2-1).getAtomBonds() + " is from " + m1bonds.get(i));
			ArrayList<Integer> m2bonds = m2parents.get(atomNum2-1).getAtomBonds();
			for (Integer k : m2bonds) {
			    possible.add(k);
			}
			for (int y = 0; y < possible.size(); y ++) {
			    int j = possible.get(y);
			    if (m2eqMatch.contains(j-1)) {
				//System.out.println(j-1 + " " + m2eqMatch);
				possible.remove(y);
				y --;
			    }
			}
			System.out.println(possible.size() + " " + possible);
			if (possible.size() == 1) {
			    int matchNum = possible.get(0);
			    ArrayList<Integer> matchBonds = m2parents.get(matchNum-1).getAtomBonds();
			    
			    //System.out.println("match");
			    mapping.put(atomNum, possible.get(0));
			    m1eqMatch.add(atomNum-1);
			    m2eqMatch.add(possible.get(0)-1);
			}
			if (possible.size() > 1) {
			    for (int q = 0; q < possible.size(); q++) {
				//System.out.println(reactantTest + " " + m2parents.get(possible.get(q)-1) + " " + possible.get(q));
				if (reactantTest.toString().equals(m2parents.get(possible.get(q)-1).toString())) {
				    mapping.put(atomNum, possible.get(q));
				    m1eqMatch.add(atomNum - 1);
				    m2eqMatch.add(possible.get(q)-1);
				}
			    }
		    }
		    }
		}
	    }
	}
	//System.out.println();
	//System.out.println();
	for (int m1atom : mapping.keySet()) {
	    System.out.println("reactant atom " + m1atom + " maps to " + mapping.get(m1atom) + " in the product");
	}
	System.out.println("Start removal testing");
	//if (addedNodes.size() == 0 && rmvedNodes.size() == 0) {
	    //addedNodes.add("C1");
	    //rmvedNodes.add("C1");
	//}
	
	// compare with removal
	for (int m1lvl1 = 0; m1lvl1 < m1parents.size(); m1lvl1++) {
	    m2loop:
	    for (int m2lvl1 = 0; m2lvl1 < m2parents.size(); m2lvl1++) {
		matchTreeNode m1root = m1parents.get(m1lvl1); // root of reactant tree for atom of index m1lvl1
		matchTreeNode m2root = m2parents.get(m2lvl1); // root of product tree for atom of index m2lvl1
		if (m1root.toString().equals(m2root.toString()) && !m1eqMatch.contains(m1lvl1) && !m2eqMatch.contains(m2lvl1)) { // root compare
		    ArrayList<matchTreeNode> m1children = m1root.getChildren(); // lvl2 of the reactant tree
		    ArrayList<matchTreeNode> m2children = m2root.getChildren(); // lvl2 of the product tree
		    ArrayList<matchTreeNode> childrenTest = compareChildren(m1children, m2children); // comparison of lvl2 of trees
		    if (childrenTest.size() == 0) { // lvl2 matches
			ArrayList<Integer> m2childrenMatched = new ArrayList<Integer>(); // indices of product lvl2 nodes matched
			for (int m1lvl2 = 0; m1lvl2 < m1children.size(); m1lvl2++) {
			    for (int m2lvl2 = 0; m2lvl2 < m2children.size(); m2lvl2++) {
				matchTreeNode m1lvl2curr = m1children.get(m1lvl2); // lvl2 node for lvl 3/4 testing (reactant)
				matchTreeNode m2lvl2curr = m2children.get(m2lvl2); // lvl2 node for lvl 3/4 testing (product)
				if (m1lvl2curr.toString().equals(m2lvl2curr.toString()) && !m2childrenMatched.contains(m2lvl2)) {
				    ArrayList<matchTreeNode> m1gchildren = m1lvl2curr.getChildren(); // reactant lvl3 nodes
				    ArrayList<matchTreeNode> m2gchildren = m2lvl2curr.getChildren(); // product lvl3 nodes
				    ArrayList<matchTreeNode> gchildrenTest = compareChildren(m1gchildren, m2gchildren); // lvl3 node comparison
				    if (gchildrenTest.size() == 0) { // lvl3 nodes match
					if (m1gchildren.size() == 0 && m2gchildren.size() == 0) {
					    m2childrenMatched.add(m2lvl2); // no gchildren and child string match so m2lvl2 matches m1lvl2.
					    if (m2childrenMatched.size() == m2children.size()) { // all m2children are matched
						int i = m1lvl1 + 1;
						int i2 = m2lvl1 + 1;
						mapping.put(i, i2);
						m1eqMatch.add(m1lvl1);
						m2eqMatch.add(m2lvl1);
						break m2loop; // m2 matched, jump back for next reactant atom matching
					    }
					} // end no gchildren block
					ArrayList<Integer> m2gchildrenMatched = new ArrayList<Integer>(); // indices of lvl3 matched via lvl4
					for (int m1lvl3 = 0; m1lvl3 < m1gchildren.size(); m1lvl3 ++) {
					    for (int m2lvl3 = 0; m2lvl3 < m2gchildren.size(); m2lvl3++) {
						matchTreeNode m1lvl3curr = m1gchildren.get(m1lvl3); // lvl3 node for testing (reactant)
						matchTreeNode m2lvl3curr = m2gchildren.get(m2lvl3); // lvl3 node for testing (product)
						if (m1lvl3curr.toString().equals(m2lvl3curr.toString()) && !m2gchildrenMatched.contains(m2lvl3)) {
						    ArrayList<matchTreeNode> m1ggchildren = m1lvl3curr.getChildren(); // reactant lvl4
						    ArrayList<matchTreeNode> m2ggchildren = m2lvl3curr.getChildren(); // product lvl4
						    ArrayList<matchTreeNode> ggchildrenTest = compareChildren(m1ggchildren, m2ggchildren); // lvl4 test
						    if (ggchildrenTest.size() == 0) { // lvl4 matches / all lvls match
							m2childrenMatched.add(m2lvl3); // lvl3 match added
							if (m2gchildrenMatched.size() == m2gchildren.size()) { // all lvl3 product matched
							    m2childrenMatched.add(m2lvl2); // lvl3 matched so lvl2 matches
							    if (m2childrenMatched.size() == m2children.size()) { // all lvl2 product matches
								int index = m1lvl1 + 1;
								int index2 = m2lvl1 + 1;
								mapping.put(index, index2);
								m1eqMatch.add(m1lvl1);
								m2eqMatch.add(m2lvl1);
								break m2loop; // test next reactant atom
							    }
							}
						    } else { // lvl4 doesnt match need removal testing
							ArrayList<matchTreeNode> m1ggmatched = new ArrayList<matchTreeNode>();
							ArrayList<matchTreeNode> ggmatched = new ArrayList<matchTreeNode>(); // product lvl4 matched 
							for (int m1lvl4 = 0; m1lvl4 < m1ggchildren.size(); m1lvl4++) {
							    m2lvl4loop:
							    for (int m2lvl4 = 0; m2lvl4 < m2ggchildren.size(); m2lvl4++) {
								matchTreeNode m1lvl4curr = m1ggchildren.get(m1lvl4); 
								matchTreeNode m2lvl4curr = m2ggchildren.get(m2lvl4);
								if (m1lvl4curr.toString().equals(m2lvl4curr.toString()) && !ggmatched.contains(m2lvl4curr)) {
								    m1ggmatched.add(m1lvl4curr);
								    ggmatched.add(m2lvl4curr);
								    break m2lvl4loop;
								}
							    } // end m2lvl4 loop for lvl4 removal
							} // end m1lvl4 loop for lvl4 removal
							if (ggmatched.size() == m2ggchildren.size()) { // double check for all matched automatically
							    m2gchildrenMatched.add(m2lvl3); // lvl4 matches so lvl3 match added
							    if (m2gchildrenMatched.size() == m2gchildren.size()) { // checks if all lvl3 are matched
								m2childrenMatched.add(m2lvl2); // lvl3 matches so add lvl2 match
								if (m2childrenMatched.size() == m2children.size()) { // lvl2 all match, add match
								    int index = m1lvl1 + 1;
								    int index2 = m2lvl1 + 1;
								    mapping.put(index, index2);
								    m1eqMatch.add(m1lvl1);
								    m2eqMatch.add(m2lvl1);
								    break m2loop;
								}
							    }
							} else { // lvl4 doesnt match, need to try removing
							    ArrayList<matchTreeNode> m1ggremoved = new ArrayList<matchTreeNode>(); // buffer for removed nodes R
							    ArrayList<matchTreeNode> m2ggremoved = new ArrayList<matchTreeNode>(); // buffer for removed nodes P
							    for (int rmv = 0; rmv < m1ggchildren.size(); rmv++) { // removes nodes removed from Reactant
								matchTreeNode m = m1ggchildren.get(rmv); // node to test
								if (rmvedNodes.contains(m.toString()) && !m1ggmatched.contains(m)) {
								    m1ggchildren.remove(m); // remove node from reactant tree
								    m1ggremoved.add(m); // add node to buffer
								}
							    }
							    for (int add = 0; add < m2ggchildren.size(); add++) { // remove nodes added to Product
								matchTreeNode m = m2ggchildren.get(add);
								if (addedNodes.contains(m.toString()) && !ggmatched.contains(m)) {
								    m2ggchildren.remove(m); // remove node from product tree
								    m2ggremoved.add(m); // add node to buffer
								}
							    }
							    ArrayList<matchTreeNode> ggretest = compareChildren(m1ggchildren, m2ggchildren); // retest ggchildren 
							    int m1atomnum = m1lvl1 + 1;
							    int m2atomnum = m2lvl1 + 1;
							    System.out.println(m1atomnum + " " + m2atomnum + " removed " + m1ggremoved + "; " + m2ggremoved + " from lvl4");
							    for (matchTreeNode mtn : m1ggremoved) {
								m1ggchildren.add(mtn); // re-add node for non-problematic testing later (reactant)
							    }
							    for (matchTreeNode mtn : m2ggremoved) {
								m2ggchildren.add(mtn); // re-add node for testing later (product)
							    }
							    if (ggretest.size() == 0) { // removal successful for lvl4
								m2gchildrenMatched.add(m2lvl3); // lvl3 match added
								if (m2gchildrenMatched.size() == m2gchildren.size()) { // all lvl3 matches
								    m2childrenMatched.add(m2lvl2); // lvl2 match added
								    if (m2childrenMatched.size() == m2children.size()) { // all lvl2 matches
									int index = m1lvl1 + 1;
									int index2 = m2lvl1 + 1;
									mapping.put(index, index2); // match
									m1eqMatch.add(m1lvl1);
									m2eqMatch.add(m2lvl1);
									break m2loop; // move to next match test
								    }
								}
							    } else { // removal in lvl4 doesnt work, so no match
								break m2loop;
							    }
							}
						    } // end gg children removal test
						} // end lvl3 String test block
						
					    } // end m2lvl3 loop
					} // end m1lvl3 loop 
				    } else { // gchildren dont match
					ArrayList<matchTreeNode> m1gchildrenMatched = new ArrayList<matchTreeNode>(); // m1gchildren which dont need removal
					ArrayList<matchTreeNode> m2gchildrenMatched = new ArrayList<matchTreeNode>(); // m2gchildren which dont need removal
					for (int m1lvl3rmv = 0; m1lvl3rmv < m1gchildren.size(); m1lvl3rmv++) {
					    for (int m2lvl3rmv = 0; m2lvl3rmv < m2gchildren.size(); m2lvl3rmv++) {
						matchTreeNode m1gchild = m1gchildren.get(m1lvl3rmv);
						matchTreeNode m2gchild = m2gchildren.get(m2lvl3rmv);
						if (m1gchild.toString().equals(m2gchild.toString()) && !m2gchildrenMatched.contains(m2gchild)) { // gchild roots match
						    ArrayList<matchTreeNode> m1ggchildren = m1gchild.getChildren();
						    ArrayList<matchTreeNode> m2ggchildren = m2gchild.getChildren();
						    ArrayList<matchTreeNode> ggchildrenTest = compareChildren(m1ggchildren, m2ggchildren);
						    if (ggchildrenTest.size() == 0) { // lvl4/3 match
							m1gchildrenMatched.add(m1gchild);
							m2gchildrenMatched.add(m2gchild);
						    } else { // lvl4 doesnt match need removal testing
							ArrayList<matchTreeNode> m1ggmatched = new ArrayList<matchTreeNode>();
							ArrayList<matchTreeNode> ggmatched = new ArrayList<matchTreeNode>(); // product lvl4 matched 
							for (int m1lvl4 = 0; m1lvl4 < m1ggchildren.size(); m1lvl4++) {
							    m2lvl4loop:
							    for (int m2lvl4 = 0; m2lvl4 < m2ggchildren.size(); m2lvl4++) {
								matchTreeNode m1lvl4curr = m1ggchildren.get(m1lvl4); 
								matchTreeNode m2lvl4curr = m2ggchildren.get(m2lvl4);
								if (m1lvl4curr.toString().equals(m2lvl4curr.toString()) && !ggmatched.contains(m2lvl4curr)) {
								    m1ggmatched.add(m1lvl4curr);
								    ggmatched.add(m2lvl4curr);
								    break m2lvl4loop;
								}
							    } // end m2lvl4 loop for lvl4 removal
							} // end m1lvl4 loop for lvl4 removal
							if (ggmatched.size() == m2ggchildren.size()) { // double check for all matched automatically
							    m2gchildrenMatched.add(m2gchild); // lvl4 matches so lvl3 match added
							    if (m2gchildrenMatched.size() == m2gchildren.size()) { // checks if all lvl3 are matched
								m2childrenMatched.add(m2lvl2); // lvl3 matches so add lvl2 match
								if (m2childrenMatched.size() == m2children.size()) { // lvl2 all match, add match
								    int index = m1lvl1 + 1;
								    int index2 = m2lvl1 + 1;
								    mapping.put(index, index2);
								    m1eqMatch.add(m1lvl1);
								    m2eqMatch.add(m2lvl1);
								    break m2loop;
								}
							    }
							} else { // lvl4 doesnt match, need to try removing
							    ArrayList<matchTreeNode> m1ggremoved = new ArrayList<matchTreeNode>(); // buffer for removed nodes R
							    ArrayList<matchTreeNode> m2ggremoved = new ArrayList<matchTreeNode>(); // buffer for removed nodes P
							    for (int rmv = 0; rmv < m1ggchildren.size(); rmv++) { // removes nodes removed from Reactant
								matchTreeNode m = m1ggchildren.get(rmv); // node to test
								if (rmvedNodes.contains(m.toString()) && !m1ggmatched.contains(m)) {
								    m1ggchildren.remove(m); // remove node from reactant tree
								    m1ggremoved.add(m); // add node to buffer
								}
							    }
							    for (int add = 0; add < m2ggchildren.size(); add++) { // remove nodes added to Product
								matchTreeNode m = m2ggchildren.get(add);
								if (addedNodes.contains(m.toString()) && !ggmatched.contains(m)) {
								    m2ggchildren.remove(m); // remove node from product tree
								    m2ggremoved.add(m); // add node to buffer
								}
							    }
							    ArrayList<matchTreeNode> ggretest = compareChildren(m1ggchildren, m2ggchildren); // retest ggchildren 
							    int m1atomnum = m1lvl1 + 1;
							    int m2atomnum = m2lvl1 + 1;
							    System.out.println(m1atomnum + " " + m2atomnum + " removed " + m1ggremoved + "; " + m2ggremoved + " from lvl4");
							    for (matchTreeNode mtn : m1ggremoved) {
								m1ggchildren.add(mtn); // re-add node for non-problematic testing later (reactant)
							    }
							    for (matchTreeNode mtn : m2ggremoved) {
								m2ggchildren.add(mtn); // re-add node for testing later (product)
							    }
							    if (ggretest.size() == 0) { // removal successful for lvl4
								m1gchildrenMatched.add(m1gchild);
								m2gchildrenMatched.add(m2gchild);
							    }
							}
						    } // end gg children removal test
						}
					    }
					}
					ArrayList<matchTreeNode> m1lvl3rmved = new ArrayList<matchTreeNode>();
					ArrayList<matchTreeNode> m2lvl3rmved = new ArrayList<matchTreeNode>();
					for (matchTreeNode m1node : m1gchildren) {
					    if (!m1gchildrenMatched.contains(m1node) && rmvedNodes.contains(m1node)) {
						m1lvl3rmved.add(m1node);
						m1gchildren.remove(m1node);
					    }
					}
					for (matchTreeNode m2node : m2gchildren) {
					    if (!m2gchildrenMatched.contains(m2node) && addedNodes.contains(m2node)) {
						m2lvl3rmved.add(m2node);
						m2gchildren.remove(m2node);
					    }
					}
					ArrayList<matchTreeNode> gchildRetest = compareChildren(m1gchildren, m2gchildren);
					int m1atomnum = m1lvl1 + 1;
					int m2atomnum = m2lvl1 + 1;
					System.out.println(m1atomnum + " " + m2atomnum + " removed " + m1lvl3rmved + "; " + m2lvl3rmved + " from lvl 3");
					for (matchTreeNode fix : m1lvl3rmved) {
					    m1gchildren.add(fix);
					}
					for (matchTreeNode fix : m2lvl3rmved) {
					    m2gchildren.add(fix);
					} 
					if (gchildRetest.size() == 0) {
					    m2childrenMatched.add(m2lvl2);
					    if (m2childrenMatched.size() == m2children.size()) {
						int index = m1lvl1 + 1;
						int index2 = m2lvl1 + 1;
						mapping.put(index, index2);
						m1eqMatch.add(m1lvl1);
						m2eqMatch.add(m2lvl1);
						break m2loop;
					    }
					}
				    } // close lvl3 removal
				} // end lvl2 node string test
			    } // end m2lvl2 loop
			} // end m1lvl2 loop			
		    } else { // lvl2 doesnt match
			ArrayList<matchTreeNode> m1rmvMatched = new ArrayList<matchTreeNode>();
			ArrayList<matchTreeNode> m2rmvMatched = new ArrayList<matchTreeNode>();
			for (int m1rmv = 0; m1rmv < m1children.size(); m1rmv++) {
			    rmvloop:
			    for (int m2rmv = 0; m2rmv < m2children.size(); m2rmv++) {
				matchTreeNode m1child = m1children.get(m1rmv);
				matchTreeNode m2child = m2children.get(m2rmv);
				if (m1child.toString().equals(m2child.toString()) && !m2rmvMatched.contains(m2child)) {
				    ArrayList<matchTreeNode> m1gchildren = m1child.getChildren();
				    ArrayList<matchTreeNode> m2gchildren = m2child.getChildren();
				    ArrayList<matchTreeNode> gchildTest = compareChildren(m1gchildren, m2gchildren);
				    if (gchildTest.size() == 0) {
					ArrayList<Integer> gchildrenMatched = new ArrayList<Integer>();
					for (int m1grmv = 0; m1grmv < m1gchildren.size(); m1grmv++) {
					    for (int m2grmv = 0; m2grmv < m2gchildren.size(); m2grmv++) {
						matchTreeNode m1gchild = m1gchildren.get(m1grmv);
						matchTreeNode m2gchild = m2gchildren.get(m2grmv);
						if (m1gchild.toString().equals(m2gchild.toString()) && !gchildrenMatched.contains(m2grmv)) {
						    ArrayList<matchTreeNode> m1ggchildren = m1gchild.getChildren();
						    ArrayList<matchTreeNode> m2ggchildren = m2gchild.getChildren();
						    ArrayList<matchTreeNode> ggtest = compareChildren(m1ggchildren, m2ggchildren);
						    if (ggtest.size() == 0) {
							gchildrenMatched.add(m2grmv);
							if (gchildrenMatched.size() == m2gchildren.size()) {
							    m1rmvMatched.add(m1child);
							    m2rmvMatched.add(m2child);
							}
						    } else { // lvl4 removal test
							ArrayList<matchTreeNode> m1ggmatched = new ArrayList<matchTreeNode>();
							ArrayList<matchTreeNode> ggmatched = new ArrayList<matchTreeNode>(); // product lvl4 matched 
							for (int m1lvl4 = 0; m1lvl4 < m1ggchildren.size(); m1lvl4++) {
							    m2lvl4loop:
							    for (int m2lvl4 = 0; m2lvl4 < m2ggchildren.size(); m2lvl4++) {
								matchTreeNode m1lvl4curr = m1ggchildren.get(m1lvl4); 
								matchTreeNode m2lvl4curr = m2ggchildren.get(m2lvl4);
								if (m1lvl4curr.toString().equals(m2lvl4curr.toString()) && !ggmatched.contains(m2lvl4curr)) {
								    m1ggmatched.add(m1lvl4curr);
								    ggmatched.add(m2lvl4curr);
								    break m2lvl4loop;
								}
							    } // end m2lvl4 loop for lvl4 removal
							} // end m1lvl4 loop for lvl4 removal
							if (ggmatched.size() == m2ggchildren.size()) { // double check for all matched automatically
							    m1rmvMatched.add(m1child);
							    m2rmvMatched.add(m2child);
							} else { // lvl4 doesnt match, need to try removing
							    ArrayList<matchTreeNode> m1ggremoved = new ArrayList<matchTreeNode>(); // buffer for removed nodes R
							    ArrayList<matchTreeNode> m2ggremoved = new ArrayList<matchTreeNode>(); // buffer for removed nodes P
							    for (int rmv = 0; rmv < m1ggchildren.size(); rmv++) { // removes nodes removed from Reactant
								matchTreeNode m = m1ggchildren.get(rmv); // node to test
								if (rmvedNodes.contains(m.toString()) && !m1ggmatched.contains(m)) {
								    m1ggchildren.remove(m); // remove node from reactant tree
								    m1ggremoved.add(m); // add node to buffer
								}
							    }
							    for (int add = 0; add < m2ggchildren.size(); add++) { // remove nodes added to Product
								matchTreeNode m = m2ggchildren.get(add);
								if (addedNodes.contains(m.toString()) && !ggmatched.contains(m)) {
								    m2ggchildren.remove(m); // remove node from product tree
								    m2ggremoved.add(m); // add node to buffer
								}
							    }
							    ArrayList<matchTreeNode> ggretest = compareChildren(m1ggchildren, m2ggchildren); // retest ggchildren 
							    int m1atomnum = m1lvl1 + 1;
							    int m2atomnum = m2lvl1 + 1;
							    System.out.println(m1atomnum + " " + m2atomnum + " removed " + m1ggremoved + "; " + m2ggremoved + " from lvl4");
							    for (matchTreeNode mtn : m1ggremoved) {
								m1ggchildren.add(mtn); // re-add node for non-problematic testing later (reactant)
							    }
							    for (matchTreeNode mtn : m2ggremoved) {
								m2ggchildren.add(mtn); // re-add node for testing later (product)
							    }
							    if (ggretest.size() == 0) { // removal successful for lvl4
								m1rmvMatched.add(m1child);
								m2rmvMatched.add(m2child);
							    }
							}
						    } // end gg children removal test
						}
						
					    }
					}
				    } else { // gchildren dont match
					ArrayList<matchTreeNode> m1gchildrenMatched = new ArrayList<matchTreeNode>(); // m1gchildren which dont need removal
					ArrayList<matchTreeNode> m2gchildrenMatched = new ArrayList<matchTreeNode>(); // m2gchildren which dont need removal
					for (int m1lvl3rmv = 0; m1lvl3rmv < m1gchildren.size(); m1lvl3rmv++) {
					    for (int m2lvl3rmv = 0; m2lvl3rmv < m2gchildren.size(); m2lvl3rmv++) {
						matchTreeNode m1gchild = m1gchildren.get(m1lvl3rmv);
						matchTreeNode m2gchild = m2gchildren.get(m2lvl3rmv);
						if (m1gchild.toString().equals(m2gchild.toString()) && !m2gchildrenMatched.contains(m2gchild)) { // gchild roots match
						    ArrayList<matchTreeNode> m1ggchildren = m1gchild.getChildren();
						    ArrayList<matchTreeNode> m2ggchildren = m2gchild.getChildren();
						    ArrayList<matchTreeNode> ggchildrenTest = compareChildren(m1ggchildren, m2ggchildren);
						    if (ggchildrenTest.size() == 0) { // lvl4/3 match
							m1gchildrenMatched.add(m1gchild);
							m2gchildrenMatched.add(m2gchild);
						    } else { // lvl4 doesnt match need removal testing
							ArrayList<matchTreeNode> m1ggmatched = new ArrayList<matchTreeNode>();
							ArrayList<matchTreeNode> ggmatched = new ArrayList<matchTreeNode>(); // product lvl4 matched 
							for (int m1lvl4 = 0; m1lvl4 < m1ggchildren.size(); m1lvl4++) {
							    m2lvl4loop:
							    for (int m2lvl4 = 0; m2lvl4 < m2ggchildren.size(); m2lvl4++) {
								matchTreeNode m1lvl4curr = m1ggchildren.get(m1lvl4); 
								matchTreeNode m2lvl4curr = m2ggchildren.get(m2lvl4);
								if (m1lvl4curr.toString().equals(m2lvl4curr.toString()) && !ggmatched.contains(m2lvl4curr)) {
								    m1ggmatched.add(m1lvl4curr);
								    ggmatched.add(m2lvl4curr);
								    break m2lvl4loop;
								}
							    } // end m2lvl4 loop for lvl4 removal
							} // end m1lvl4 loop for lvl4 removal
							if (ggmatched.size() == m2ggchildren.size()) { // double check for all matched automatically
							    m1rmvMatched.add(m1child);
							    m2rmvMatched.add(m2child);
							} else { // lvl4 doesnt match, need to try removing
							    ArrayList<matchTreeNode> m1ggremoved = new ArrayList<matchTreeNode>(); // buffer for removed nodes R
							    ArrayList<matchTreeNode> m2ggremoved = new ArrayList<matchTreeNode>(); // buffer for removed nodes P
							    for (int rmv = 0; rmv < m1ggchildren.size(); rmv++) { // removes nodes removed from Reactant
								matchTreeNode m = m1ggchildren.get(rmv); // node to test
								if (rmvedNodes.contains(m.toString()) && !m1ggmatched.contains(m)) {
								    m1ggchildren.remove(m); // remove node from reactant tree
								    m1ggremoved.add(m); // add node to buffer
								}
							    }
							    for (int add = 0; add < m2ggchildren.size(); add++) { // remove nodes added to Product
								matchTreeNode m = m2ggchildren.get(add);
								if (addedNodes.contains(m.toString()) && !ggmatched.contains(m)) {
								    m2ggchildren.remove(m); // remove node from product tree
								    m2ggremoved.add(m); // add node to buffer
								}
							    }
							    ArrayList<matchTreeNode> ggretest = compareChildren(m1ggchildren, m2ggchildren); // retest ggchildren 
							    int m1atomnum = m1lvl1 + 1;
							    int m2atomnum = m2lvl1 + 1;
							    System.out.println(m1atomnum + " " + m2atomnum + " removed " + m1ggremoved + "; " + m2ggremoved + " from lvl4");
							    for (matchTreeNode mtn : m1ggremoved) {
								m1ggchildren.add(mtn); // re-add node for non-problematic testing later (reactant)
							    }
							    for (matchTreeNode mtn : m2ggremoved) {
								m2ggchildren.add(mtn); // re-add node for testing later (product)
							    }
							    if (ggretest.size() == 0) { // removal successful for lvl4
								m1gchildrenMatched.add(m1gchild);
								m2gchildrenMatched.add(m2gchild);
							    }
							}
						    } // end gg children removal test
						}
					    }
					}
					ArrayList<matchTreeNode> m1lvl3rmved = new ArrayList<matchTreeNode>();
					ArrayList<matchTreeNode> m2lvl3rmved = new ArrayList<matchTreeNode>();
					for (matchTreeNode m1node : m1gchildren) {
					    if (!m1gchildrenMatched.contains(m1node) && rmvedNodes.contains(m1node)) {
						m1lvl3rmved.add(m1node);
						m1gchildren.remove(m1node);
					    }
					}
					for (matchTreeNode m2node : m2gchildren) {
					    if (!m2gchildrenMatched.contains(m2node) && addedNodes.contains(m2node)) {
						m2lvl3rmved.add(m2node);
						m2gchildren.remove(m2node);
					    }
					}
					ArrayList<matchTreeNode> gchildRetest = compareChildren(m1gchildren, m2gchildren);
					int m1atomnum = m1lvl1 + 1;
					int m2atomnum = m2lvl1 + 1;
					System.out.println(m1atomnum + " " + m2atomnum + " removed " + m1lvl3rmved + "; " + m2lvl3rmved + " from lvl 3");
					for (matchTreeNode fix : m1lvl3rmved) {
					    m1gchildren.add(fix);
					}
					for (matchTreeNode fix : m2lvl3rmved) {
					    m2gchildren.add(fix);
					} 
					if (gchildRetest.size() == 0) {
					    m1rmvMatched.add(m1child);
					    m2rmvMatched.add(m2child); 
					}
				    } // close lvl3 removal
				    
				}
			    }
			} // found children for removal testing
			ArrayList<matchTreeNode> m1rmvedChildren = new ArrayList<matchTreeNode>();
			ArrayList<matchTreeNode> m2rmvedChildren = new ArrayList<matchTreeNode>();
			for (matchTreeNode m1chldrmv : m1children) {
			    if (rmvedNodes.contains(m1chldrmv) && !m1rmvMatched.contains(m1chldrmv)) {
				m1rmvedChildren.add(m1chldrmv);
				m1children.remove(m1chldrmv);
			    }
			}
			for (matchTreeNode m2chldrmv : m2children) {
			    if (addedNodes.contains(m2chldrmv) && !m2rmvMatched.contains(m2chldrmv)) {
				m2rmvedChildren.add(m2chldrmv);
				m2children.remove(m2chldrmv);
			    }
			}
			ArrayList<matchTreeNode> chldrentest = compareChildren(m1children, m2children);
			int m1atomnum = m1lvl1 + 1;
			int m2atomnum = m2lvl1 + 1;
			System.out.println(m1atomnum + " " + m2atomnum + " removed " + m1rmvedChildren + "; " + m2rmvedChildren + " from lvl2");
			for (matchTreeNode n : m1rmvedChildren) {
			    m1children.add(n);
			}
			for (matchTreeNode n : m2rmvedChildren) {
			    m2children.add(n);
			}
			if (chldrentest.size() == 0) {
			    int index = m1lvl1 + 1;
			    int index2 = m2lvl1 + 1;
			    mapping.put(index, index2);
			    m1eqMatch.add(m1lvl1);
			    m2eqMatch.add(m2lvl1);
			}
		    } // end lvl2 remove test
		} // end m1root and m2root string test
	    } // end m2lvl1 loop
	} // end m1lvl1 loop
	

	// post removal double checking for left over nodes if no bond changes
	if (addedNodes.size() == 0 && rmvedNodes.size() == 0) {
	    for (int m1parIndex = 0; m1parIndex < m1parents.size(); m1parIndex++) {
		matchTreeNode reactantTest = m1parents.get(m1parIndex);
		ArrayList<Integer> m1bonds = reactantTest.getAtomBonds();
		int atomNum = reactantTest.getAtomNumber();
		//System.out.println("Reactant atom " + atomNum + " is bonded to " + m1bonds);
		ArrayList<Integer> possible = new ArrayList<Integer>();
		for (int i=0; i<m1bonds.size(); i++){
		    if (mapping.keySet().contains(m1bonds.get(i)) && !mapping.keySet().contains(atomNum)) {
			int bondedAtom = mapping.get(m1bonds.get(i));
			int atomNum2 = m2parents.get(bondedAtom).getAtomNumber()-1;
			//System.out.println("Product atom " + atomNum2 + " is bonded to " + m2parents.get(atomNum2-1).getAtomBonds() + " is from " + m1bonds.get(i));
			ArrayList<Integer> m2bonds = m2parents.get(atomNum2-1).getAtomBonds();
			for (Integer k : m2bonds) {
			    possible.add(k);
			}
			for (int y = 0; y < possible.size(); y ++) {
			    int j = possible.get(y);
			    if (m2eqMatch.contains(j-1)) {
				//System.out.println(j-1 + " " + m2eqMatch);
				possible.remove(y);
				y --;
			    }
			}
			if (possible.size() == 1) {
			    //System.out.println("match");
			    mapping.put(atomNum, possible.get(0));
			    m1eqMatch.add(atomNum-1);
			    m2eqMatch.add(possible.get(0)-1);
			} if (possible.size() > 1) {
			    for (int q = 0; q < possible.size(); q++) {
				//System.out.println(reactantTest + " " + m2parents.get(possible.get(q)-1) + " " + possible.get(q));
				if (reactantTest.toString().equals(m2parents.get(possible.get(q)-1).toString())) {
				    mapping.put(atomNum, possible.get(q));
				    m1eqMatch.add(atomNum - 1);
				    m2eqMatch.add(possible.get(q)-1);
				}
			    }
			}
		    }
		}
	    }
	}
    	System.out.println(mapping.keySet().size());
	for (int m1atom : mapping.keySet()) {
	    System.out.println("reactant atom " + m1atom + " maps to " + mapping.get(m1atom) + " in the product");
	}
	
    } // close matchAtoms method

   

    public ArrayList<matchTreeNode> compareChildren(ArrayList<matchTreeNode> m1, ArrayList<matchTreeNode> m2) {
	ArrayList<matchTreeNode> unMatchedNodes = new ArrayList<matchTreeNode>();
	ArrayList<Integer> m2visited = new ArrayList<Integer>();
	if (m1.size() != m2.size()) {
	    unMatchedNodes.add(new matchTreeNode("z",4,null,0, new ArrayList<String>()));
	    return unMatchedNodes;
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
	for (int k=0; k<m2.size(); k++) {
	    if (!m2visited.contains(k)) {
		unMatchedNodes.add(m2.get(k));
	    }
	}
	return unMatchedNodes;
    } // close compareChildren method
}