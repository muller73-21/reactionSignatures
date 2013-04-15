import java.util.*;
import java.io.*;

public class MoleculeCreationTest {
    public static void main (String args[]) throws FileNotFoundException {
	ArrayList<Molecule> reactants = new ArrayList<Molecule>();
	ArrayList<Molecule> products = new ArrayList<Molecule>();
	MoleculeCreationTest test = new MoleculeCreationTest();
	MoleculeCreationTest test2 = new MoleculeCreationTest();
	Scanner reactantSc = new Scanner(new File(args[0]));
	Molecule reactant = test.createMolecule (reactantSc);
	reactants.add(reactant);
	System.out.println(reactants.get(0));
	System.out.println(reactants.get(0).getBonds());
	System.out.println(reactant.getBonds());
	
	Scanner productSc = new Scanner(new File(args[1]));
	System.out.println("scanner creation " + reactant.getBonds());
	Molecule product = test2.createMolecule(productSc);
	System.out.println(" molecule creation " + reactant.getBonds());
	products.add(product);
	//reactants.add(product);
	
	System.out.println(reactants.get(0));
	System.out.println(products.get(0));

	System.out.println("Reactant = " + reactants.get(0).getBonds());
	System.out.println(reactant.getBonds());
	//System.out.println("Product = " + reactants.get(1).getBonds());
	System.out.println("Product  = " + products.get(0).getBonds());
	System.out.println(product.getBonds());
    }
    public Molecule createMolecule(Scanner molsc){
	System.out.println("reached createMolecule");
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
}