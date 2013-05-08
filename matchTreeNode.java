import java.util.*;
import java.io.*;

/***
    Class is created to be able to build atom match trees, for purpose of atom
    matching between reactants and products.
 */
public class matchTreeNode {
    private String atomSym;
    private int bondType;
    private ArrayList<matchTreeNode> children;
    private matchTreeNode parent;
    private int atomNumber;
    
    /***
	Constructor without children
	@param s string representing atom symbol
	@param i number of bonds between parent atom and passed atom
	@param p node which is the parent of this node in the tree
	@param num the number of the atom in the mol file
     */
    public matchTreeNode(String s, int i, matchTreeNode p, int num) {
	atomSym = s;
	bondType = i;
	parent = p;
	children = new ArrayList<matchTreeNode>();
	atomNumber = num;
    }

    /***
	Constructor with children
	@param s string representing atom symbol
	@param i number of bonds between parent atom and passed atom
	@param p node which is the parent of this node in the tree
	@param list current children to be the children of the node being created.
	@param num number of atom in the mol file.
     */
    public matchTreeNode(String s, int i, matchTreeNode p, ArrayList<matchTreeNode> list, int num) {
	atomSym = s;
	bondType = i;
	parent = p;
	children = list;
	atomNumber = num;
    }

    /***
	@return returns number of atom in the mol file 
     */
    public int getAtomNumber() {
	return atomNumber;
    }
    
    /***
	@return String with format of bondtype directly following atom symbol in a string
     */
    public String toString() {
	return atomSym + bondType;
    }

    /***
	@return string which is the symbol of the atom in this node
     */
    public String getAtom() {
	return atomSym;
    }

    /***
	@return bondtype number 1,2,3,etc depending on single, double, triple bond, etc.
     */
    public int getBondType() {
	return bondType;
    }
    
    /***
	@return list of children nodes for the given node.
     */
    public ArrayList<matchTreeNode> getChildren() {
	return children;
    }

    /***
	@param chldrn list of new children nodes
     */
    public void setChildren(ArrayList<matchTreeNode> chldrn) {
	children = chldrn;
    }
    
    /***
	@return the tree node which is the parent of this current node
     */
    public matchTreeNode getParent() {
	return parent;
    }
    
    /***
	@param p node to set to parent in case where one wants to change parent of node
     */
    public void setParent(matchTreeNode p) {
	parent = p;
    }

    /***
	@return true if node has no children (aka is a leaf of tree), false if it has children
     */
    public boolean isLeaf() {
	if (children.size() == 0) 
	    return true;
	return false;
    }

    /***
	@return true if node has no parent, aka is root of tree, false if parent exists.
     */
    public boolean isRoot() {
	if (parent == null)
	    return true;
	return false;
    }

    /***
	Adds new child to node's child list, sets child's parent to this node automatically
	@param child tree node which wants to be added as a child to this node
     */
    public void addChild(matchTreeNode child) {
	child.setParent(this);
	children.add(child);
    }

    /***
	Adds child with just passing atomSym and bondType, automatically sets parent to this node.
	@param atom symbol for atom which is child (bonded) to this atom
	@param bond type of bond, single, double, etc.
     */
    public void addChild(String atom, int bond, int atomNum) {
	children.add(new matchTreeNode(atom, bond, this, atomNum));
    }
}