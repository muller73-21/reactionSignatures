import java.io.*;
import java.util.*;

/***
    Representation of a Molecule including atom list, 
number of bonds, number of atoms and a list of string representations
of bonds within the molecule. 
 */

public class Molecule {
    private ArrayList<String> atoms = new ArrayList<String>();
    private int numberOfBonds;
    private int numberOfAtoms;
    private ArrayList<String> bonds = new ArrayList<String>();
    private int[] atomBondFreqs;
    private int numOfMolecules;
    private ArrayList<String> charges = new ArrayList<String>();
    private HashMap<String, Integer> bondTypes = new HashMap<String, Integer>();
    private ArrayList<ArrayList<String>> matchTrees = new ArrayList<ArrayList<String>>();
    private ArrayList<matchTreeNode> atomTrees = new ArrayList<matchTreeNode>();
    
    /***
	Constructs a new Molecule Object - molecule objects can hold information
	about more than 1 molecule, given that all molecules were contained within
	1 mol file. Constructs empty atomBondFreqs array with size numberOfAtoms.
	
	@param atms List of atoms contained within this molecule
	@param bonds number of bonds in the molecule (double and triple bonds count only as one).
	@param numAtms number of atoms within the molecule, should be equal to size of List atms passed in
	@param bnds List of String representation of bonds in molecule with the structure "x y z" where x 
	and y are number of the atom in the list of atoms (1 indexed) and z is the type of bond. x,y,z are integers.
    */
    public Molecule(ArrayList<String> atms, int bonds, int numAtms, 
		    ArrayList<String> bnds, ArrayList<String> chrgs) {
	atoms = atms;
	numberOfBonds = bonds;
        numberOfAtoms = numAtms;
	this.bonds = bnds;
	atomBondFreqs = new int[numberOfAtoms];
	charges = chrgs;
    }
    /***
	Method populates the freq array with the number of bonds each atom has. Does this by
	scanning the bonds list and for each x and y take the corresponding spot in the array 
	and increases the number by z. This will be useful to tell how many hydrogens are attached to each
	atom in the molecule if the hydrogens are not explicitly in the mol file. Generates seperate from
	the time of Molecule construction so that the frequencies can be recalculated if the atoms are rearranged
	and the bond numbers change.
     */
    public void generateFreq() {
	for (String bond : bonds) {
	    Scanner sc = new Scanner(bond);
	    int firstAtom = Integer.parseInt(sc.next());
	    int secondAtom = Integer.parseInt(sc.next());
	    int freq = Integer.parseInt(sc.next());
	    atomBondFreqs[firstAtom - 1] += freq;
	    atomBondFreqs[secondAtom - 1] += freq;
	}
    }
    /***
	Method returns the atomBondFreqs array.
	@return atomBondFreqs number of atoms on each molecule, index of atom matches index in atom List..
     */
    public int[] getFreqs() {
	return atomBondFreqs;
    }
    /***
	return number of bonds on a certain atom.
	@param index of atom in the atom List
	@return returns freq number at given index.
     */
    public int getFreqOfAtom(int index) {
	return atomBondFreqs[index];
    }
    /***
	Replaces the bonds List in the Molecule with new bonds List
	@param bnds new list of bonds "x y z" which will replace the bonds of the molecule.
     */
    public void setBonds(ArrayList<String> bnds) {
	bonds = bnds;
    }
    /***
	Gets the List of bonds of the molecule
	@return List of bonds in the molecule of format "x y z".
     */
    public ArrayList<String> getBonds() {
	return bonds;
    }
    /***
	sets the number of bonds to the inputted bonds
	@param bonds number of bonds to replace the old number of bonds.
     */
    public void setNumOfBonds(int bonds) {
	numberOfBonds = bonds;
    }
    /***
	Return numberOfBonds in the Molecule
	@return number of Bonds in the Molecule.
     */
    public int getNumberOfBonds () {
	return numberOfBonds;
    }
    /***
	Replace the numberOfAtoms in Molecule with new number
	@param atoms number of Atoms want to be in the molecule.
     */
    public void setNumOfAtoms(int atoms) {
	numberOfAtoms = atoms;
    }
    /***
	Return number of Atoms in the molecule.
	@return numberOfAtoms in molecule.
     */
    public int getNumberOfAtoms() {
	return numberOfAtoms;
    }
    /***
	Replaces Atom List with new Atom List
	@param atms List of Atoms to replace current Atoms in Molecule Object
     */
    public void setAtoms(ArrayList<String> atms){
	atoms = atms;
    }
    /***
	Return List of Match Trees for atoms in molecule
	@return List of Match Trees for atoms in molecule
     */
    public ArrayList<ArrayList<String>> getMatchTrees() {
	return matchTrees;
    }
    /***
	Return List of Atoms in Molecule
	@return List of Atoms in Molecule Object.
     */
    public ArrayList<String> getAtoms() {
	return atoms;
    }
    /***
	Return list containing parent nodes of match trees
	@return List of matchTreeNodes which are heads of atom Match trees
     */
    public ArrayList<matchTreeNode> getatomTrees() {
	return atomTrees;
    }

    /***
	Use connectivity trees to detect multiple molecules in mol file
	@return number of molecules in mol file which created molecule object
     */
    public int findNumOfMlcs() {
	int[] connections = new int[numberOfAtoms];
	for (int i = 0; i<connections.length; i++) {
	    connections[i] = -1;
	}
	for (String s : bonds) {
	    Scanner sc = new Scanner(s);
	    int a = Integer.parseInt(sc.next());
	    int b = Integer.parseInt(sc.next());
	    if (a < b) {
		connections[b-1] = a-1;
	    } else {
		connections[a-1] = b-1;
	    } 
	}
	numOfMolecules = 0;
	for (int i=0; i<connections.length; i++) {
	    if (connections[i] == -1) {
		numOfMolecules ++;
	    }
	}
	return numOfMolecules;
    }
    
    /***
	Creates a map of the types of bonds in molecule
	mapped to frequency at which they occur. 
	@return mapping of bond types to frequency of bond type.
     */
    public HashMap<String, Integer> listChangedBonds() {
	for (String s : bonds) {
	    Scanner sc = new Scanner(s);
	    String atomA = atoms.get(Integer.parseInt(sc.next()) - 1);
	    String atomB = atoms.get(Integer.parseInt(sc.next()) -1);
	    int bondtype = Integer.parseInt(sc.next());
	    String bond = "";
	    String rbond = "";
	    if (bondtype == 1) {
		bond = bond + atomA + " - " + atomB;
		rbond = rbond + atomB + " - " + atomA;
	    } else if (bondtype == 2) {
		bond = bond + atomA + " = " + atomB;
		rbond = rbond + atomB + " = " + atomA;
	    } else if (bondtype == 3) {
		bond = bond + atomA + " t " + atomB;
		rbond = rbond + atomB + " t " + atomA;
	    }
	    if (bondTypes.containsKey(bond)) {
		int curFreq = bondTypes.get(bond);
		curFreq++;
		bondTypes.remove(bond);
		bondTypes.put(bond, curFreq);
	    } else if (bondTypes.containsKey(rbond)) {
		int curFreq = bondTypes.get(rbond);
		curFreq++;
		bondTypes.remove(rbond);
		bondTypes.put(rbond, curFreq);
	    } else {
		bondTypes.put(bond, 1);
	    }
	}
	return bondTypes;
    }
    
    public void buildMatchTrees() {
	matchTreeNode atom = new matchTreeNode(atoms.get(0), 0, null, 1, this.bonds);
	ArrayList<matchTreeNode> parentNodes = new ArrayList<matchTreeNode>();
	
	for (int i = 0; i<atoms.size(); i++) {
	    matchTrees.add(new ArrayList<String>());
	    int currentAtom = i+1;
	    atom = new matchTreeNode(atoms.get(i), 0, null, i+1, this.bonds);
	    parentNodes.add(atom);
	    for (String bond : bonds) {
		Scanner sc = new Scanner(bond);
		//System.out.print(bond);
		int atom1 = Integer.parseInt(sc.next());
		int atom2 = Integer.parseInt(sc.next());
		int bondtype = Integer.parseInt(sc.next());
		//System.out.println(atoms.get(atom1-1) + " " + atoms.get(atom2-1));
		if (atom1 == currentAtom || atom2 == currentAtom) {
		    if (atom1 == currentAtom) {
			atom.addChild(atoms.get(atom2 - 1), bondtype, atom2);
		    } else {
			atom.addChild(atoms.get(atom1 - 1), bondtype, atom1);
		    }
		}
	    }
	    for (matchTreeNode child : atom.getChildren()) {
		buildChildren(child);
		for (matchTreeNode grandchild : child.getChildren()) {
		    buildChildren(grandchild);
		}
	    }
	    
	    int index = i + 1;
	    String root = atom.toString();
	    String path = "";
	    atomTrees = parentNodes;
	    //System.out.println(atom + " " + index );
	    for (matchTreeNode child : atom.getChildren()) {
		//System.out.print(child + ": ");
		path = root + child.toString();
		if (child.getChildren().size() == 0) 
		    matchTrees.get(i).add(path);
		for (matchTreeNode grandchild : child.getChildren()) {
		    //System.out.print(grandchild + " {");
		    path = path + grandchild;
		    if (grandchild.getChildren().size() == 0) 
			matchTrees.get(i).add(path);
		    for (matchTreeNode ggrandchild : grandchild.getChildren()) {
			//System.out.print(ggrandchild + " ");
			matchTrees.get(i).add(path + ggrandchild);
		    }
		    //System.out.print("} ");
		    path = root + child.toString();
		}
		//System.out.println();
		path = root;
	    }
	}
    }

    private void buildChildren(matchTreeNode atom) {
	int currentAtom = atom.getAtomNumber();
	int parentAtom = atom.getParent().getAtomNumber();
	for (String bond : bonds) {
	    Scanner sc = new Scanner(bond);
	    //System.out.println(bond);
	    int atom1 = Integer.parseInt(sc.next());
	    int atom2 = Integer.parseInt(sc.next());
	    int bondtype = Integer.parseInt(sc.next());
	    //System.out.println(atoms.get(atom1-1) + " " + atoms.get(atom2-1));
	    if (atom1 == currentAtom && atom2 != parentAtom) {
		atom.addChild(new matchTreeNode(atoms.get(atom2 - 1), bondtype, atom, atom2, this.bonds));
	    } else if (atom2 == currentAtom && atom1 != parentAtom) {
		atom.addChild(new matchTreeNode(atoms.get(atom1 - 1), bondtype, atom, atom1, this.bonds));
	    }
	    
	}
    }
}


