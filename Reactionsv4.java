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
	if (addedNodes.size() == 0 && rmvedNodes.size() == 0) {
	    //addedNodes.add("C1");
	    //rmvedNodes.add("C1");
	}
	
	// compare with removal
	
	for (int m1lvl1 = 0; m1lvl1 < m1parents.size(); m1lvl1 ++) {
	    lvl1continue:
	    for (int m2lvl1 = 0; m2lvl1 < m2parents.size(); m2lvl1++) {
		matchTreeNode m1root = m1parents.get(m1lvl1);
		matchTreeNode m2root = m2parents.get(m2lvl1);
		if (m1root.toString().equals(m2root.toString()) && !m1eqMatch.contains(m1lvl1) && !m2eqMatch.contains(m2lvl1)) {
		    ArrayList<matchTreeNode> m1children = m1root.getChildren();
		    ArrayList<matchTreeNode> m2children = m2root.getChildren();
		    ArrayList<matchTreeNode> childMatch = compareChildren(m1children, m2children);
		    
		    if (childMatch.size() == 0 ) { // children match
			
		    } else { // children don't match
			ArrayList<Integer> m1chMatched = new ArrayList<Integer>();
			ArrayList<Integer> m2chMatched = new ArrayList<Integer>();
			for (int m1lvl2=0; m1lvl2<m1children.size(); m1lvl2++) {
			    lvl2continue:
			    for (int m2lvl2=0; m2lvl2<m2children.size(); m2lvl2++) {
				matchTreeNode m1child = m1children.get(m1lvl2);
				matchTreeNode m2child = m2children.get(m2lvl2);
				if (m1child.toString().equals(m2child.toString()) && !m2chMatched.contains(m2lvl2)) {
				    ArrayList<matchTreeNode> m1gchildren = m1child.getChildren();
				    ArrayList<matchTreeNode> m2gchildren = m2child.getChildren();
				    ArrayList<matchTreeNode> gchildMatch = compareChildren(m1gchildren, m2gchildren);
				    if (gchildMatch.size() == 0) {
					ArrayList<Integer> m2gchMatched = new ArrayList<Integer>();
					for (int m1lvl3=0; m1lvl3<m1gchildren.size(); m1lvl3++) {
					    lvl3continue:
					    for (int m2lvl3=0; m2lvl3<m2gchildren.size(); m2lvl3++) {
						matchTreeNode m1gchild = m1gchildren.get(m1lvl3);
						matchTreeNode m2gchild = m2gchildren.get(m2lvl3);
						if (m1gchild.toString().equals(m2gchild.toString()) && !m2gchMatched.contains(m2lvl3)) {
						    ArrayList<matchTreeNode> m1ggchildren = m1gchild.getChildren();
						    ArrayList<matchTreeNode> m2ggchildren = m2gchild.getChildren();
						    ArrayList<matchTreeNode> ggchildMatch = compareChildren(m1ggchildren, m2ggchildren);
						    if (ggchildMatch.size() == 0) {
							m2gchMatched.add(m2lvl3);
							if (m2gchMatched.size() == m2gchildren.size()) {
							    m2chMatched.add(m2lvl2);
							    m1chMatched.add(m1lvl2);
							}
						    } else {
							break lvl3continue;
						    }
						} else {
						    break lvl3continue;
						}
					    }
					}
				    } else {
					break lvl2continue;
				    }
				} else {
				    break lvl2continue;
				}
			    }
			} // found children who don't match
			ArrayList<matchTreeNode> m1childrenRemoved = new ArrayList<matchTreeNode>();
			ArrayList<matchTreeNode> m2childrenRemoved = new ArrayList<matchTreeNode>();
			for (int m1rem = 0 ; m1rem < m1children.size(); m1rem++) {
			    matchTreeNode m1child = m1children.get(m1rem);
			    if (!m1chMatched.contains(m1rem) && rmvdNodes.contains(m1child.toString())) {
				m1children.remove(m1child);
				m1childrenRemoved.add(m1child);
			    }
			}
			for (int m2rem = 0; m2rem < m2children.size(); m2rem++) {
			    matchTreeNode m2child = m2children.get(m2rem);
			    if (!m2chMatched.contains(m2rem) && addedNodes.contains(m2child.toString())) {
				m2children.remove(m2child);
				m2childrenRemoved.add(m2child);
			    }
			}
			if (m1chMatched.size() + m1childrenRemoved.size() == m1children.size()) {
			    if (m2chMatched.size() + m2childrenRemoved.size() == m2children.size()) {
				ArrayList<matchTreeNode> childrenRetest = compareChildren(m1children, m2children);
				if (childrenRetest.size() == 0) {
				    
				}
			    }
			}
		    } // end children removal
		} else {
		    break lvl1continue;
		}
	    } // close m2lvl2 loop
	} // close m1lvl1 loop

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