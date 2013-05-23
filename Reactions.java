import java.io.*;
import java.util.*;
/***
 *  Main Program running to generate Reaction Signatures from given mol files
 *  Currently takes mol files as a commandline parameter in the format:
 *  {reactants} break {products} - the word break separating reactant from product 
 *  must be spelled out in all lower case letters. 
 *  @author Matt Muller
 */
public class Reactions {
    /***
	ArrayList of type Molecule holding reactant molecule objects built from 
	reactant mol files (those before "break"). 
     */
    public ArrayList<Molecule> reactants = new ArrayList<Molecule>();
    /***
	ArrayList of type Molecule holding product molecule objects built from 
	product mol files (those after "break").
     */
    public ArrayList<Molecule> products = new ArrayList<Molecule>();
    /***
	Main method, passes command line parameters to createReaction() to 
	fill the ArrayLists with the appropriate molecule objects.
	Cycles through reactant Molecules and generates frequencies of atoms 
	in bonds to judge how many bonds explicitly said in mol files are attached
	to each atom. Finds number of molecules in each mol file, compares the bond of 
	each mol file to find out how many bonds change and of what type they are. 

	@param list of mol files to be processed with the word "break" in lower case letters 
	separating the reactant mol files from the product mol files.
     */
    public static void main (String args[]) throws FileNotFoundException {
	Reactions reactions = new Reactions();
	reactions.createReaction(args);
	for (Molecule mlc : reactions.reactants) {
	    for (int i = 0; i< mlc.getAtoms().size(); i++) {
		System.out.print(mlc.getAtoms().get(i) + " ");
	    }
	    System.out.println();
	    mlc.generateFreq();
	    int [] freq = mlc.getFreqs();
	    for (int i=0; i<freq.length; i++) {
		System.out.print(freq[i] + " ");
	    }
	    System.out.println();
	    System.out.println("molecules = " + mlc.findNumOfMlcs());
	    HashMap<String, Integer> bondTypes = mlc.listChangedBonds();
	    Set<String> bondPics = bondTypes.keySet();
	    for (String s :bondPics) {
		System.out.println(s + " " + bondTypes.get(s));
	    }
	    mlc.buildMatchTrees();
	    int atomNumCount = 1;
	    for (ArrayList<String> atomTree : mlc.getMatchTrees()) {
		System.out.print(atomNumCount + ": ");
		for (String s : atomTree) {
		    System.out.print(s + ", ");
		} 
		System.out.println();
		atomNumCount ++;
	    }
	} 
	for (Molecule mlc : reactions.products) {
	    mlc.generateFreq();
	    int [] freq = mlc.getFreqs();
	    for (int i=0; i<freq.length; i++) {
		System.out.print(freq[i] + " ");
	    }
	    System.out.println();
	    System.out.println("molecules = " + mlc.findNumOfMlcs());
	    HashMap<String, Integer> bondTypes = mlc.listChangedBonds();
	    Set<String> bondPics = bondTypes.keySet();
	    for (String s :bondPics) {
		System.out.println(s + " " + bondTypes.get(s));
	    }
	    mlc.buildMatchTrees();
	    int atomNumCount = 1;
	    for (ArrayList<String> atomTree : mlc.getMatchTrees()) {
		System.out.print(atomNumCount + ": ");
		for (String s : atomTree) {
		    System.out.print(s + ", ");
		} 
		System.out.println();
		atomNumCount ++;
	    }
	}
	reactions.matchAtoms(reactions.reactants.get(0), reactions.products.get(0));
    }
    
    /***
	Method takes a list of filenames with the word "break" seperating reactants from
	products. Loops through the list creating Molecule object for each mol file before 
	the word "break" and populates the ArrayList reactants with them in order. Catches 
	the "break" and then creates Molecule objects for each mol file after, populating 
	those into the ArrayList products. Afterwards, prints the Bonds, Atoms, NumberOfAtoms,
	and NumOfBonds fields of the Molecules in each list. This is for purpose of testing that
	the molecule object creation matches the mol files passed in.

	@param args list of mol file names with reactants and products separated with "break".
     */
    public void createReaction(String[] args) throws FileNotFoundException {
	File reactant;
	Scanner reactantSc;
	Molecule molecule;
	File product;
	Scanner productSc;
	Boolean isAReactantFile = true;
	for (int i=0; i<args.length; i++) {
	    if (args[i].equals("break")) {
		isAReactantFile = false;
		i++;
	    }
	    if (isAReactantFile == true) {
		reactant = new File(args[i]);
		reactantSc = new Scanner (reactant);
		molecule  = createMolecule(reactantSc);
		reactants.add(molecule);
	    }
	    if (isAReactantFile == false){
		System.out.println(args[i]);
		product = new File(args[i]);
		productSc = new Scanner(product);
		Molecule rmolecule  = createMolecule(productSc);
		products.add(rmolecule);
	    }	    
	}
	int i = 0;
	for (Molecule m : reactants) {
	    i++;
	    System.out.println("Reactant " + i + " # of atoms = " + m.getNumberOfAtoms());
	    System.out.println("Reactant " + i + " # of bonds = " + m.getNumberOfBonds());
	    System.out.println("Reactant " + i + " atoms = " + m.getAtoms());
	    System.out.println("Reactant " + i + " bonds = " + m.getBonds());
	}
	i = 0;
	for (Molecule m : products) {
	    i++;
	    System.out.println("Product " + i + " # of atoms = " + m.getNumberOfAtoms());
	    System.out.println("Product " + i + " # of bonds = " + m.getNumberOfBonds());
	    System.out.println("Product " + i + " atoms = " + m.getAtoms());
	    System.out.println("Product " + i + " bonds = " + m.getBonds());
	}
    }
    /***
	Method takes a scanner object built around a molfile and creates a Molecule object
	using data extracted from the molfile. This allows for the easy access to specific 
	information contained within the mol file without needing constant processing of the
	file itself. 

	@param molsc Scanner built around a file of type .mol 
	@return Molecule object built from information contained in the scanner's mol file.
     */
    public Molecule createMolecule(Scanner molsc){
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
	while (temp.charAt(2) == ' '){
	    Scanner tempsc = new Scanner(temp);
	    String atom = tempsc.next();
	    atom = tempsc.next();
	    atom = tempsc.next();
	    atom = tempsc.next();
	    atoms.add(atom);
	    temp = molsc.nextLine();
	}
	Scanner bondsc = new Scanner(temp);
	while (bondsc.hasNextInt()){
	    String bond = "" + bondsc.nextInt() + " " + bondsc.nextInt() + " " + bondsc.nextInt();	    
	    bonds.add(bond);
	    temp = molsc.nextLine();
	    bondsc = new Scanner (temp);
	}
	Molecule molecule = new Molecule(atoms, numberOfBonds, numberOfAtoms, bonds);
	return molecule;
    }

    /***

     */
    public void matchAtoms(Molecule m1, Molecule m2) {
	Map<Integer, Integer> mapping = new HashMap<Integer, Integer>();
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
		    }else {
			int key1Freq = m1BondTypes.get(key);
			int key2Freq = m2BondTypes.get(key);
			if (key1Freq - key2Freq > 0) {
			    changedBonds.get("-").add(key);
			} else if (key1Freq - key2Freq < 0) {
			    changedBonds.get("+").add(key);
			}
		    }
		}
		firstThrough = true;
	    }
	}
	ArrayList<String> chngdNodes = new ArrayList<String>();
	ArrayList<String> addedNodes = new ArrayList<String>();
	for (String add : changedBonds.get("+")) {
	    System.out.println("+ " + add);
	    String parentAtom = "" + add.charAt(0);
	    String bondedAtom = "" + add.charAt(2);
	    String bondtype = "" + add.charAt(1);
	    if (bondtype.equals("-")) {
		
		addedNodes.add( bondedAtom + "1");
	
	    } else if (bondtype.equals("=")) {
	
		addedNodes.add(bondedAtom + "2");
	
	    } else {
		
		addedNodes.add( bondedAtom + "3");
	
	    }
	}
	for (String lose : changedBonds.get("-")) {
	    System.out.println("- " + lose);
	    String parentAtom = "" + lose.charAt(0);
	    String bondedAtom = "" + lose.charAt(2);
	    String bondtype = "" + lose.charAt(1);
	    if (bondtype.equals("-")) {
		chngdNodes.add( bondedAtom + "1");
	
	    } else if (bondtype.equals("=")) {
		chngdNodes.add( parentAtom + "2");
	
	    } else {
		chngdNodes.add( bondedAtom + "3");
		
	    }
	}
	System.out.println("Added nodes " + addedNodes);
	ArrayList<ArrayList<String>> m1atomTrees = m1.getMatchTrees();
	ArrayList<ArrayList<String>> m2atomTrees = m2.getMatchTrees();
	ArrayList<Integer> m2AtomsMapped = new ArrayList<Integer>();
	ArrayList<Integer> m2PathsUsed = new ArrayList<Integer>();
	ArrayList<Integer> m1PathsMapped = new ArrayList<Integer>();
	HashMap<Integer, Integer> atomMap = new HashMap<Integer, Integer>();
	ArrayList<matchTreeNode> m1parents = m1.getatomTrees();
	ArrayList<matchTreeNode> m2parents = m2.getatomTrees();
	ArrayList<Integer> m2visited = new ArrayList<Integer>();
	for (int m1parIndex = 0; m1parIndex < m1parents.size(); m1parIndex++) {
	    matchTreeNode m1curr = m1parents.get(m1parIndex);
	    for (int m2parIndex = 0; m2parIndex < m2parents.size(); m2parIndex++) {
		matchTreeNode m2curr = m2parents.get(m2parIndex);
		if (m1curr.toString().equals(m2curr.toString()) && !m2visited.contains(m2parIndex)) {
		    ArrayList<matchTreeNode> m1children1 = m1curr.getChildren();
		    ArrayList<matchTreeNode> m2children1 = m2curr.getChildren();
		    ArrayList<matchTreeNode> match = compareChildren(m1children1, m2children1);
		    //System.out.println("Match " + match + " size " + match.size());
		    if (match.size() != 0) {		    
			int changeMatch = 0;
			for (matchTreeNode m : match) {
			    //System.out.println(m + " " + addedNodes.contains(m.toString()));
			    if (addedNodes.contains(m.toString())){
				
				changeMatch++;
			    }
			}
			//System.out.println(changeMatch);
			if (changeMatch == match.size()) {
			    for (matchTreeNode m: match) {
				m2children1.remove(m);
				//System.out.println("removed something");
			    }
			}
			m2curr.setChildren(m2children1);
			//System.out.println("chgndNodes = " + chngdNodes);
			matchTreeNode diff = null;
			ArrayList<matchTreeNode> toBeRemoved = new ArrayList<matchTreeNode>();
			for (matchTreeNode n : m1children1) {
			    //System.out.println(n);
			    if ( chngdNodes.contains(n.toString())) {
				//System.out.println(n);
				//m1children1.remove(n);
				toBeRemoved.add(n);
				//System.out.println(m1children1);
			    }
			}
		       
			for (int remove = 0; remove < toBeRemoved.size();remove++) {
			    m1children1.remove(toBeRemoved.get(remove));
			}
			m1curr.setChildren(m1children1);
			//System.out.println(m1children1 + "; " + m2curr.getChildren());
			match = compareChildren(m1children1,m2curr.getChildren());
			if (match.size() != 0 || m1children1.size() != m2children1.size()) {
			    continue;
			} else {
			    int index = m1parIndex + 1;
			    int ind = m2parIndex + 1;
			    //System.out.println(m1curr.toString() +" " + index  +  ": "+ m1children1 + " "  + m2curr.toString() + " " + ind + ": " + m2children1 + " " + match);
			    mapping.put(index, ind);
			    m2visited.add(m2parIndex);
			    break;
			}
			
		    } else {
			int index = m1parIndex + 1;
			int ind = m2parIndex + 1;
			//System.out.println(m1curr.toString() +" " + index +  ": "+ m1children1 + " "  + m2curr.toString() + " " + ind + ": " + m2children1 + " " + match);
			mapping.put(index, ind);
			m2visited.add(m2parIndex);
			break;
		    }
		    
		} 
	    }
	}
	/*for (int m1Index = 0; m1Index < m1atomTrees.size(); m1Index++) { // cycles through each atom in m1
	    atomMap = new HashMap<Integer, Integer>();
	    for (int m2Index = 0; m2Index < m2atomTrees.size(); m2Index ++) { // cycles through each atom in m2
		if (!m2AtomsMapped.contains(m2Index + 1)) {
		    if (m1atomTrees.get(m1Index).size() == m2atomTrees.get(m2Index).size()) {			
			for (int m1pathIndex = 0; m1pathIndex < m1atomTrees.get(m1Index).size(); m1pathIndex++) { //go through each path in atom selected from outer Loop in m1			    
			    m2PathsUsed = new ArrayList<Integer>();
			    for (int m2pathIndex = 0; m2pathIndex < m2atomTrees.get(m2Index).size(); m2pathIndex++) { // goes through each tree in m2
				if (m1atomTrees.get(m1Index).get(m1pathIndex).equals(m2atomTrees.get(m2Index).get(m2pathIndex)) && !atomMap.containsValue(m2pathIndex)) {
				    atomMap.put(m1pathIndex, m2pathIndex);
				    m2PathsUsed.add(m2pathIndex);
				    break;
				}
			    }			    
			    if (atomMap.keySet().size() == m1atomTrees.get(m1Index).size()){
				mapping.put(m1Index + 1, m2Index + 1);
				m2AtomsMapped.add(m2Index + 1);
			    }
			}
		    }
		    if (m2AtomsMapped.contains(m2Index+1)) {
			break;
		    }		   
		}		
	    }
	}
	*/
	System.out.println(mapping.keySet().size());
	for (int m1atom : mapping.keySet()) {
	    System.out.println("reactant atom " + m1atom + " maps to " + mapping.get(m1atom)+ " in the product");
	}
    }
    
    private ArrayList<matchTreeNode> compareChildren(ArrayList<matchTreeNode> m1, ArrayList<matchTreeNode> m2) {
	ArrayList<matchTreeNode> unMatchedNodes = new ArrayList<matchTreeNode>();
	ArrayList<Integer> m2visited = new ArrayList<Integer>();
	/*if (m1.size() != m2.size()) {
	   return false;
	   }*/
	
	for (int i=0; i<m1.size(); i++) {	   
	    String temp = m1.get(i).toString();
	    boolean match = false;
	    for (int j=0; j<m2.size();j++) {
		if (temp.equals(m2.get(j).toString()) && !m2visited.contains(j)) {
		    m2visited.add(j);
		    match = true;
		    break;
		} else {
		    continue;
		}
	    } 
	    if (match == false) {
		//unMatchedNodes.add(m1.get(i));
		//return false;
	    }	    
	}
	for (int k=0;k<m2.size();k++) {
	    if (!m2visited.contains(k)) {
		unMatchedNodes.add(m2.get(k));
	    }
	}
	return unMatchedNodes;
	
    }
}