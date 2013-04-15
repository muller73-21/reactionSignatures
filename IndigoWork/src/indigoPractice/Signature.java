package indigoPractice;
import java.util.*;
import java.io.*;

import com.ggasoftware.indigo.Indigo;
import com.ggasoftware.indigo.IndigoObject;

import com.ggasoftware.indigo.*;


@SuppressWarnings("unused")
public class Signature {
	
	
	public static String molFileName;
	public static Indigo indigo;
	
	public static void main(String[] args) throws FileNotFoundException  {
		File f = new File("2-propone.mol");
		indigo = new Indigo();
		//IndigoObject mol1 = indigo.loadMolecule("C1CCC(=O)CC1");
		//IndigoObject mol2 = indigo.loadMolecule("C1CCC(CC1)O");
		IndigoObject mol2 = indigo.loadMoleculeFromFile("2-propone.mol");
		/*System.out.println("Atoms in molecule = " + mol1.countAtoms());
		System.out.println("Bonds in molecule = " + mol1.countBonds());
		System.out.println("Molecule is: " + mol1.grossFormula());
		System.out.println(mol1.molfile());
		*/
		System.out.println("Atoms in molecule = " + mol2.countAtoms());
		System.out.println("Bonds in molecule = " + mol2.countBonds());
		System.out.println("Molecule is: " + mol2.grossFormula());
		System.out.println(mol2.molfile());
		/*
		File molFile1 = new File("phenol.txt");
		System.out.println(molFile1);
		Scanner sc = new Scanner(molFile1);
		while (sc.hasNext()){
			System.out.println(sc.next());
		}
		*/
	}
}
