/***
    Version 5 of reactions taking version 4 revising the removal comparisons and trying to add in working prematch comparisons for all tests.
*/
import java.io.*;
import java.util.*;

public class Reactionsv5 {
    /***
	Molecule representing mol file of all reactants for reaction
     */
    public Molecule reactants;
    /***
	Molecule representing mol file of all products for reaction
     */
    public Molecule products;
    /***
	Main method creating reaction and running matching algorithm and eventually signature creation algorithm.
    */
    public static void main(String[] args) throws FileNotFoundException {
	Reactionsv5 reactions = new Reactionsv5();
	reactions.createReaction(args);
	Molecule reactant = reactions.reactants;
	for (int i = 0; i < reactant.getAtoms().size(); i++) {
	    System.out.print(reactant.getAtoms().get(i) + " ");
	}
	System.out.println();
	reactant.generateFreq();
	int [] freq = reactant.getFreqs();
	for (int i=0; i<freq.length; i++) {
	    System.out.print(freq[i] + " " );
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
	    System.out.print(atomNumCount + ": " );
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
	    System.out.print(patomNumCount + ": " );
	    for (String s : patomTree) {
		System.out.print(s + ", ");
	    }
	    System.out.println();
	    patomNumCount++;
	}
	reactions.matchAtoms(reactant, product);
    }
    /***
	Method which takes mol files and converts them into molecule objects 
     */
    public void createReaction(String[] args) throws FileNotFoundException {
	File reactant;
	Scanner reactantSc;
	Molecule rmolecule;
	File product;
	Scanner productSc;
	Molecule pmolecule;
	if (args.length != 2) {
	    System.out.println("Wrong number of arguments, enter 2 arguments: reactant product");
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
    /***
	method taking scanner and creating molecule object to create reaction representation
     */
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
	Scanner infosc = new Scanner(info);
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
	while(bondsc.hasNextInt()) {
	    String bond = "" + bondsc.nextInt() + " " + bondsc.nextInt() + " " + bondsc.nextInt();
	    bonds.add(bond);
	    temp = molsc.nextLine();
	    bondsc = new Scanner(temp);
	}
	Molecule molecule = new Molecule(atoms, numberOfBonds, numberOfAtoms, bonds);
	return molecule;
    }
    
    public void matchAtoms (Molecule m1, Molecule m2) {
	Map<Integer, Integer> mapping = new HashMap <Integer, Integer> ();
	int m1AtomNum = 0;
	int m2AtomNum = 0;
	int rootMatches = 0;
	HashMap<String, Integer> m1BondTypes = m1.listChangedBonds();
	HashMap<String, Integer> m2BondTypes = m2.listChangedBonds();
	Set<String> m1keys = m1BondTypes.keySet();
	HashMap<String, ArrayList<String>> changedBonds = new HashMap<String, ArrayList<String>>();
	changedBonds.put("+", new ArrayList<String>());
	changedBonds.put("-", new ArrayList<String>());
	boolean firstThrough = false;
	for (String key : m1keys) {
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
		}// close key 2 for loop
	    }
	} // close key for loop
	
	ArrayList<String> rmvedNodes = new ArrayList<String>();
	ArrayList<String> addedNodes = new ArrayList<String>();
	for (String add : changedBonds.get("+")) {
	    Scanner bondsc = new Scanner (add);
	    String parentAtom = bondsc.next();
	    String bondType = bondsc.next();
	    String bondedAtom = bondsc.next();
	    if (bondType.equals("-")) {
		addedNodes.add(bondedAtom + "1");
		addedNodes.add(parentAtom + "1");
	    } else if (bondType.equals("=")) {
		addedNodes.add(bondedAtom + "2"); // try adding reversing double bond too?? 
		addedNodes.add(parentAtom + "2");
	    } else {
		addedNodes.add(bondedAtom + "3"); // try adding reversing triple bond??
		addedNodes.add(parentAtom + "3");
	    }
	}
	for (String lose : changedBonds.get("-")) {
	    Scanner bondsc = new Scanner(lose);
	    String parentAtom = bondsc.next();
	    String bondtype = bondsc.next();
	    String bondedAtom = bondsc.next();
	    if (bondtype.equals("-")) {
		rmvedNodes.add(bondedAtom + "1"); // test added double remved nodes possibly??
		rmvedNodes.add(parentAtom + "1");
	    } else if (bondtype.equals("=")) {
		rmvedNodes.add(bondedAtom + "2"); // double adding??
		rmvedNodes.add(parentAtom + "2");
	    } else {
		rmvedNodes.add(bondedAtom + "3");
		rmvedNodes.add(bondedAtom + "3");
	    }
	}
	System.out.println("Added Nodes: " + addedNodes);
	System.out.println("Lost Nodes: " + rmvedNodes);
	HashMap<Integer, Integer> atomMap = new HashMap<Integer, Integer>();
	ArrayList<Integer> m2visited = new ArrayList<Integer>();
	ArrayList<Integer> m2eqMatch = new ArrayList<Integer>();
	ArrayList<Integer> m1eqMatch = new ArrayList<Integer>();
	ArrayList<matchTreeNode> m1parents = m1.getatomTrees();
	ArrayList<matchTreeNode> m2parents = m2.getatomTrees();
	// Start equality testing loops
	for (int m1lvl1 = 0; m1lvl1 < m1parents.size(); m1lvl1++) {
	    m2loop:
	    for (int m2lvl1 = 0; m2lvl1 < m2parents.size(); m2lvl1++) {
		matchTreeNode m1root = m1parents.get(m1lvl1);
		matchTreeNode m2root = m2parents.get(m2lvl1);
		if (m1root.toString().equals(m2root.toString()) && !m2eqMatch.contains(m2lvl1)) {
		    ArrayList<matchTreeNode> m1children = m1root.getChildren();
		    ArrayList<matchTreeNode> m2children = m2root.getChildren();
		    ArrayList<matchTreeNode> childrenTest = compareChildren(m1children, m2children);
		    if (childrenTest.size() == 0) {
			ArrayList<Integer> m2childrenMatched = new ArrayList<Integer>(); 
			for (int m1lvl2 = 0; m1lvl2 < m1children.size(); m1lvl2++) {
			    for (int m2lvl2 = 0; m2lvl2 < m2children.size(); m2lvl2++) {
				matchTreeNode m1lvl2curr = m1children.get(m1lvl2);
				matchTreeNode m2lvl2curr = m2children.get(m2lvl2);
				if (m1lvl2curr.toString().equals(m2lvl2curr.toString()) && !m2childrenMatched.contains(m2lvl2)) {
				    ArrayList<matchTreeNode> m1gchildren = m1lvl2curr.getChildren();
				    ArrayList<matchTreeNode> m2gchildren = m2lvl2curr.getChildren();
				    ArrayList<matchTreeNode> gchildrenTest = compareChildren(m1gchildren, m2gchildren);
				    if (gchildrenTest.size() == 0) {
					if (m1gchildren.size() == 0 && m2gchildren.size() == 0) {
					    m2childrenMatched.add(m2lvl2);
					    if (m2childrenMatched.size() == m2children.size()) {
						int i = m1lvl1 + 1;
						int i2 = m2lvl1 + 1;
						mapping.put(i, i2);
						m1eqMatch.add(m1lvl1);
						m2eqMatch.add(m2lvl1);
						break m2loop; // m2 matched jump back to next m1lvl1
					    }
					} // no gchildren (aka tree 2 levels deep)
					ArrayList<Integer> m2gchildrenMatched = new ArrayList<Integer>();
					for (int m1lvl3 = 0; m1lvl3 < m1gchildren.size(); m1lvl3++) {
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
								break m2loop; // break back to start of m2lvl1 loop
							    }
							}
						    }
						}
					    } // close equality m2lvl3 loop
					} // close equality m1lvl3 loop
				    }
				} 
			    } // close equality m2lvl2 loop
			} // close equality m1lvl2 loop
		    }
		}
	    } // close equality m2lvl1 loop
	} // close equality m1lvl1 loop
	System.out.println(mapping.keySet().size());
	for (int m1atom : mapping.keySet()) {
	    System.out.println("reactant atom " + m1atom + " maps to " + mapping.get(m1atom) + " in the product");
	}
	// test equality based testing??

	// test removal testing
	for (int m1lvl1 = 0; m1lvl1 < m1parents.size(); m1lvl1++) {
	    m2loop:
	    for (int m2lvl1 = 0; m2lvl1 < m2parents.size(); m2lvl1++) {
		matchTreeNode m1root = m1parents.get(m1lvl1);
		matchTreeNode m2root = m2parents.get(m2lvl1);
		if (m1root.toString().equals(m2root.toString()) && !m2eqMatch.contains(m2lvl1) && !m1eqMatch.contains(m1lvl1)) {
		    ArrayList<matchTreeNode> m1children = m1root.getChildren();
		    ArrayList<matchTreeNode> m2children = m2root.getChildren();
		    ArrayList<matchTreeNode> childrenTest = compareChildren(m1children, m2children);
		    if (childrenTest.size() == 0) {
			ArrayList<Integer> m2childrenMatched = new ArrayList<Integer>(); 
			System.out.println(m1lvl1 + " " + m2lvl1 + " children match!!!! " + m2childrenMatched.size());				
			for (int m1lvl2 = 0; m1lvl2 < m1children.size(); m1lvl2++) {
			    for (int m2lvl2 = 0; m2lvl2 < m2children.size(); m2lvl2++) {
				matchTreeNode m1lvl2curr = m1children.get(m1lvl2);
				matchTreeNode m2lvl2curr = m2children.get(m2lvl2);
				if (m1lvl2curr.toString().equals(m2lvl2curr.toString()) && !m2childrenMatched.contains(m2lvl2)) {
				    ArrayList<matchTreeNode> m1gchildren = m1lvl2curr.getChildren();
				    ArrayList<matchTreeNode> m2gchildren = m2lvl2curr.getChildren();
				    ArrayList<matchTreeNode> gchildrenTest = compareChildren(m1gchildren, m2gchildren);
				    if (gchildrenTest.size() == 0) {
					if (m1gchildren.size() == 0 && m2gchildren.size() == 0) {
					    m2childrenMatched.add(m2lvl2);
					    if (m2childrenMatched.size() == m2children.size()) {
						int i = m1lvl1 + 1;
						int i2 = m2lvl1 + 1;
						mapping.put(i, i2);
						m2eqMatch.add(m2lvl1);
						m1eqMatch.add(m1lvl1);
						break m2loop; // m2 matched jump back to next m1lvl1
					    }
					} // no gchildren (aka tree 2 levels deep)
					ArrayList<Integer> m2gchildrenMatched = new ArrayList<Integer>();
					System.out.println(m1lvl1 + " " + m2lvl1 + " Gchildren match!!!! " + m2childrenMatched.size() + " " + m2gchildrenMatched.size());

					for (int m1lvl3 = 0; m1lvl3 < m1gchildren.size(); m1lvl3++) {
					    for (int m2lvl3 = 0; m2lvl3 < m2gchildren.size(); m2lvl3++) {
						matchTreeNode m1lvl3curr = m1gchildren.get(m1lvl3);
						matchTreeNode m2lvl3curr = m2gchildren.get(m2lvl3);
						if (m1lvl3curr.toString().equals(m2lvl3curr.toString()) && !m2gchildrenMatched.contains(m2lvl3)) {
						    ArrayList<matchTreeNode> m1ggchildren = m1lvl3curr.getChildren();
						    ArrayList<matchTreeNode> m2ggchildren = m2lvl3curr.getChildren();
						    ArrayList<matchTreeNode> ggchildrenTest = compareChildren(m1ggchildren, m2ggchildren);
						    if (ggchildrenTest.size() == 0) {
							System.out.println(m1lvl1 + " " + m2lvl1 + " GGchildren match!!! " + m2childrenMatched.size() + " " + m2gchildrenMatched.size());
							m2gchildrenMatched.add(m2lvl3);
							if (m2gchildrenMatched.size() == m2gchildren.size()) {
							    m2childrenMatched.add(m2lvl2);
							    if (m2childrenMatched.size() == m2children.size()) {
								int index = m1lvl1 + 1;
								int index2 = m2lvl1 + 1;
								mapping.put(index, index2);
								m2eqMatch.add(m2lvl1);
								m1eqMatch.add(m1lvl1);
								break m2loop; // break back to start of m2lvl1 loop
							    }
							}
						    } else {
							ArrayList<matchTreeNode> m1ggmatched = new ArrayList<matchTreeNode>();
							ArrayList<matchTreeNode> ggmatched = new ArrayList<matchTreeNode>();
							for (int m1lvl4 = 0; m1lvl4 < m1ggchildren.size(); m1lvl4++) {
							    m2lvl4loop:
							    for (int m2lvl4 = 0; m2lvl4 < m2ggchildren.size(); m2lvl4++) {
								matchTreeNode m1lvl4curr = m1ggchildren.get(m1lvl4);
								matchTreeNode m2lvl4curr = m2ggchildren.get(m2lvl4);
								if (m1lvl4curr.toString().equals(m2lvl4curr.toString()) && !ggmatched.contains(m2lvl4curr)) {
								    ggmatched.add(m2lvl4curr);
								    break m2lvl4loop;
								}
							    }
							} // end finding of matched nodes
							if (ggmatched.size() == m2ggchildren.size()) {
							    m2gchildrenMatched.add(m2lvl3);
							    if (m2gchildrenMatched.size() == m2gchildren.size()) {
								m2childrenMatched.add(m2lvl2);
								if (m2childrenMatched.size() == m2children.size()) {
								    int index = m1lvl1 + 1;
								    int index2 = m2lvl1 + 1;
								    mapping.put(index, index2);
								    m2eqMatch.add(m2lvl1);
								    m1eqMatch.add(m1lvl1);
								    break m2loop;
								}
							    }
							} else {
							    ArrayList<matchTreeNode> m1ggremoved = new ArrayList<matchTreeNode>();
							    ArrayList<matchTreeNode> m2ggremoved = new ArrayList<matchTreeNode>();
							    for (int rmv = 0; rmv < m1ggchildren.size(); rmv++) {
								matchTreeNode m = m1ggchildren.get(rmv);
								if (rmvedNodes.contains(m.toString()) && !m1ggmatched.contains(m)) {
								    m1ggchildren.remove(m);
								    m1ggremoved.add(m);
								    System.out.println(m1lvl1 + " " + m2lvl1 + " removing gg " + m + " from m1");
								}
							    }
							    for (int add = 0; add < m2ggchildren.size(); add++) {
								matchTreeNode m = m2ggchildren.get(add);
								if (addedNodes.contains(m.toString()) && !ggmatched.contains(m)) {
								    m2ggchildren.remove(m);
								    m2ggremoved.add(m);
								    System.out.println(m1lvl1 + " " + m2lvl1 + " removing gg " + m + " from m2");
								}
							    }
							    ArrayList<matchTreeNode> ggretest = compareChildren(m1ggchildren, m2ggchildren);
							    if (ggretest.size() == 0) {
								System.out.println(m1lvl1 + " " + m2lvl1 + " GGretest match!!");
								m2gchildrenMatched.add(m2lvl3);
								for (matchTreeNode mtn : m1ggremoved) {
								    m1ggchildren.add(mtn);
								}
								for (matchTreeNode mtn : m2ggremoved) {
								    m2ggchildren.add(mtn);
								}
								if (m2gchildrenMatched.size() == m2gchildren.size()) {
								    m2childrenMatched.add(m2lvl2);
								    if (m2childrenMatched.size() == m2children.size()) {
									int index = m1lvl1 + 1;
									int index2 = m2lvl1 + 1;
									mapping.put(index, index2);
									m2eqMatch.add(m2lvl1);
									m1eqMatch.add(m1lvl1);
									break m2loop;
								    }
								}
							    } else {
								for (matchTreeNode mtn : m1ggremoved) {
								    m1ggchildren.add(mtn);
								}
								for (matchTreeNode mtn : m2ggremoved) {
								    m2ggchildren.add(mtn);
								}
								//break m2loop;
							    }
							}
							
						    } // end ggchildren removal
						}
					    } // close remove m2lvl3 loop
					} // close remove m1lvl3 loop
				    } // end gchildrentest == 0 loop
				    else {
					ArrayList<matchTreeNode> m1gchildrenMatch = new ArrayList<matchTreeNode>();
					ArrayList<matchTreeNode> m2gchildrenMatch = new ArrayList<matchTreeNode>();
					for (int m1gchild = 0; m1gchild < m1gchildren.size(); m1gchild++) {
					    m2gloop:
					    for (int m2gchild = 0; m2gchild < m2gchildren.size(); m2gchild++) {
						matchTreeNode m1gcurr = m1gchildren.get(m1gchild);
						matchTreeNode m2gcurr = m2gchildren.get(m2gchild);
						if (m1gcurr.toString().equals(m2gcurr.toString()) && !m2gchildrenMatch.contains(m2gcurr)) {
						    ArrayList<matchTreeNode> m1ggch = m1gcurr.getChildren();
						    ArrayList<matchTreeNode> m2ggch = m2gcurr.getChildren();
						    ArrayList<matchTreeNode> ggchtest = compareChildren(m1ggch, m2ggch);
						    if (ggchtest.size() == 0) {
							m1gchildrenMatch.add(m1gcurr);
							m2gchildrenMatch.add(m2gcurr);
							break m2gloop;
						    } else {
							break m2gloop;
						    }
						}
					    }
					}
					ArrayList<matchTreeNode> m1gremoved = new ArrayList<matchTreeNode> ();
					ArrayList<matchTreeNode> m2gremoved = new ArrayList<matchTreeNode> ();
					for (int rmv = 0; rmv < m1gchildren.size(); rmv ++) {
					    matchTreeNode m = m1gchildren.get(rmv);
					    if (rmvedNodes.contains(m.toString()) && !m1gchildrenMatch.contains(m)) {
						m1gchildren.remove(m);
						m1gremoved.add(m);
						System.out.println(m1lvl1 + " " + m2lvl1 + " removing g " + m + " from m1");
					    }
					}
					for (int add = 0; add < m2gchildren.size(); add++) {
					    matchTreeNode m = m2gchildren.get(add);
					    if (addedNodes.contains(m.toString()) && !m2gchildrenMatch.contains(m)) {
						m2gchildren.remove(m);
						m2gremoved.add(m);
						System.out.println(m1lvl1 + " " + m2lvl1 + " removing g " + m + " from m2");
					    }
					}
					ArrayList<matchTreeNode> gchildretest = compareChildren(m1gchildren, m2gchildren);
					if (gchildretest.size() == 0) {
					    System.out.println(m1lvl1 + " " + m2lvl1 + " Gchildren retest match!!!");
					    m2childrenMatched.add(m2lvl2);
					    for (matchTreeNode mtn : m1gremoved ) {
						m1gchildren.add(mtn);
					    }
					    for (matchTreeNode mtn : m2gremoved ) {
						m2gchildren.add(mtn);
					    }
					    if (m2childrenMatched.size() == m2children.size()) {
						int index1 = m1lvl1 + 1;
						int index2 = m2lvl1 + 1;
						mapping.put(index1, index2);
						m2eqMatch.add(m2lvl1);
						m1eqMatch.add(m1lvl1);
						break m2loop;
					    }
					} else {
					    for (matchTreeNode mtn : m1gremoved) {
						m1gchildren.add(mtn);
					    }
					    for (matchTreeNode mtn : m2gremoved) {
						m2gchildren.add(mtn);
					    }
					    //break m2loop;
					}
					
				    } // end gchildren removal
				} 
			    } // close removal m2lvl2 loop
			} // close removal m1lvl2 loop
		    } // end childrenTest.size() == 0 
		    else {
			ArrayList<matchTreeNode> m1childMatch = new ArrayList<matchTreeNode> ();
			ArrayList<matchTreeNode> m2childMatch = new ArrayList<matchTreeNode> ();
			for (int m1child = 0; m1child < m1children.size(); m1child++) {
			    m2childloop:
			    for (int m2child = 0; m2child < m2children.size(); m2child++) {
				matchTreeNode m1ch = m1children.get(m1child);
				matchTreeNode m2ch = m2children.get(m2child);
				if (m1ch.toString().equals(m2ch.toString()) && !m1childMatch.contains(m1ch) && !m2childMatch.contains(m2ch)) {
				    ArrayList<matchTreeNode> m1gch = m1ch.getChildren();
				    ArrayList<matchTreeNode> m2gch = m2ch.getChildren();
				    ArrayList<matchTreeNode> gchTest = compareChildren(m1gch, m2gch);
				    if (gchTest.size() == 0) {
					ArrayList<Integer> gchMatch = new ArrayList<Integer>();
					for (int m1g =0; m1g<m1gch.size(); m1g++) {
					    for (int m2g = 0; m2g < m2gch.size(); m2g++) {
						matchTreeNode m1gchcurr = m1gch.get(m1g);
						matchTreeNode m2gchcurr = m2gch.get(m2g);
						if (m1gchcurr.toString().equals(m2gchcurr.toString()) && !gchMatch.contains(m2g)) {
						    ArrayList<matchTreeNode> m1ggch = m1gchcurr.getChildren();
						    ArrayList<matchTreeNode> m2ggch = m2gchcurr.getChildren();
						    ArrayList<matchTreeNode> ggchTest = compareChildren(m1ggch, m2ggch);
						    if (ggchTest.size() == 0) {
							gchMatch.add(m2g);
							if (gchMatch.size() == m2gch.size()) {
							    m1childMatch.add(m1ch);
							    m2childMatch.add(m2ch);
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
			ArrayList<matchTreeNode> m1removed = new ArrayList<matchTreeNode>();
			ArrayList<matchTreeNode> m2removed = new ArrayList<matchTreeNode>();
			for (int i = 0; i< m1children.size(); i++) {
			    matchTreeNode m = m1children.get(i);
			    if (rmvedNodes.contains(m.toString()) && !m1childMatch.contains(m)) {
				m1children.remove(m);
				m1removed.add(m);
				System.out.println(m1lvl1 + " " + m2lvl1 + " removing c " + m + " from m1");
			    }
			}
			for (int i = 0; i<m2children.size(); i++) {
			    matchTreeNode m = m2children.get(i);
			    if (addedNodes.contains(m.toString()) && !m2childMatch.contains(m)) {
				m2children.remove(m);
				m2removed.add(m);
				System.out.println(m1lvl1 + " " + m2lvl1 + " removing c " + m + " from m2");
			    }
			}
			ArrayList<matchTreeNode> childrenreTest = compareChildren(m1children, m2children);
			if (childrenreTest.size() == 0) {
			    for (matchTreeNode mtn : m1removed ) {
				m1children.add(mtn);
			    }
			    for (matchTreeNode mtn : m2removed) {
				m2children.add(mtn);
			    }
			    System.out.println(m1lvl1 + " " + m2lvl1 + " childrenRetest match!!!");
			    int index1 = m1lvl1 + 1;
			    int index2 = m2lvl1 + 1;
			    mapping.put(index1, index2);
			    m1eqMatch.add(m1lvl1);
			    m2eqMatch.add(m2lvl1);
			    break m2loop;
			} else {
			    for (matchTreeNode mtn : m1removed) {
				m1children.add(mtn);
			    }
			    for (matchTreeNode mtn : m2removed) {
				m2children.add(mtn);
			    }
			    //break m2loop;
			}
		    }
		}
	    } // close removal m2lvl1 loop
	} // close removal m1lvl1 loop
	System.out.println(mapping.keySet().size());
	for (int m1atom : mapping.keySet()) {
	    System.out.println("reactant atom " + m1atom + " maps to " + mapping.get(m1atom) + " in the product");
	}
    }

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