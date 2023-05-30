package org.cisrg.mapping;

import java.util.ArrayList;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;


/**
 * Removal of listeners
 * 
 * @author edmund
 *
 */
public class QueryAtomContainer2 extends QueryAtomContainer {
	
	 /**
     *  Constructs an empty AtomContainer.
     */
    public QueryAtomContainer2(IChemObjectBuilder builder) {
        this(10, 10, 0, 0, builder);
    }
    
    /**
     *  Constructs an empty AtomContainer that will contain a certain number of
     *  atoms and electronContainers. It will set the starting array lengths to the
     *  defined values, but will not create any Atom or ElectronContainer's.
     *
     *@param  atomCount        Number of atoms to be in this container
     *@param  bondCount        Number of bonds to be in this container
     *@param  lpCount          Number of lone pairs to be in this container
     *@param  seCount          Number of single electrons to be in this container
     *
     */
    public QueryAtomContainer2(int atomCount, int bondCount, int lpCount, int seCount, IChemObjectBuilder builder) {
        super(builder);
        this.atomCount = 0;
        this.bondCount = 0;
        this.lonePairCount = 0;
        this.singleElectronCount = 0;
        atoms = new IAtom[atomCount];
        bonds = new IBond[bondCount];
        lonePairs = new ILonePair[lpCount];
        singleElectrons = new ISingleElectron[seCount];
        stereoElements = new ArrayList<IStereoElement>(atomCount / 2);
    }

	
	 public QueryAtomContainer2(IAtomContainer container, IChemObjectBuilder builder) {
	        super(builder);
	        this.atomCount = container.getAtomCount();
	        this.bondCount = container.getBondCount();
	        this.lonePairCount = container.getLonePairCount();
	        this.singleElectronCount = container.getSingleElectronCount();
	        this.atoms = new IAtom[this.atomCount];
	        this.bonds = new IBond[this.bondCount];
	        this.lonePairs = new ILonePair[this.lonePairCount];
	        this.singleElectrons = new ISingleElectron[this.singleElectronCount];

	        stereoElements = new ArrayList<IStereoElement>(atomCount / 2);

	        for (int f = 0; f < container.getAtomCount(); f++) {
	            atoms[f] = container.getAtom(f);
	            atoms[f].removeListener(container);	            
	            //container.getAtom(f).addListener(this);
	        }
	        for (int f = 0; f < this.bondCount; f++) {
	            bonds[f] = container.getBond(f);
	            //container.getBond(f).addListener(this);
	        }
	        for (int f = 0; f < this.lonePairCount; f++) {
	            lonePairs[f] = container.getLonePair(f);
	            //container.getLonePair(f).addListener(this);
	        }
	        for (int f = 0; f < this.singleElectronCount; f++) {
	            singleElectrons[f] = container.getSingleElectron(f);
	            //container.getSingleElectron(f).addListener(this);
	        }
	    }
	 
	 /**
	     *  Sets the array of atoms of this AtomContainer.
	     *
	     *@param  atoms  The array of atoms to be assigned to this AtomContainer
	     *@see           #getAtom
	     */
	    @Override
	    public void setAtoms(IAtom[] atoms) {
	        this.atoms = atoms;
	        //for (IAtom atom : atoms) {
	          //  atom.addListener(this);
	        //}
	        this.atomCount = atoms.length;
	        notifyChanged();
	    }
	    
	    
	    /**
	     * Sets the array of bonds of this AtomContainer.
	     *
	     * @param  bonds  The array of bonds to be assigned to
	     *                             this AtomContainer
	     * @see  #getBond
	     */
	    @Override
	    public void setBonds(IBond[] bonds) {
	        this.bonds = bonds;
	        //for (IBond bond : bonds) {
	        //    bond.addListener(this);
	        //}
	        this.bondCount = bonds.length;
	    }

	    
	    /**
	     *  Sets the atom at position <code>number</code> in [0,..].
	     *
	     *@param  number  The position of the atom to be set.
	     *@param  atom    The atom to be stored at position <code>number</code>
	     *@see            #getAtom(int)
	     */
	    @Override
	    public void setAtom(int number, IAtom atom) {
	        //atom.addListener(this);
	        atoms[number] = atom;
	        notifyChanged();
	    }

	    
	    
	    /**
	     *  Adds an atom to this container.
	     *
	     *@param  atom  The atom to be added to this container
	     */
	    @Override
	    public void addAtom(IAtom atom) {
	        if (contains(atom)) {
	            return;
	        }

	        if (atomCount + 1 >= atoms.length) {
	            growAtomArray();
	        }
	        //atom.addListener(this);
	        atoms[atomCount] = atom;
	        atomCount++;
	        notifyChanged();
	    }
	    
	    
	    /**
	     *  Grows the atom array by a given size.
	     *
	     *@see    #growArraySize
	     */
	    private void growAtomArray() {
	        growArraySize = (atoms.length < growArraySize) ? growArraySize : atoms.length;
	        IAtom[] newatoms = new IAtom[atoms.length + growArraySize];
	        System.arraycopy(atoms, 0, newatoms, 0, atoms.length);
	        atoms = newatoms;
	    }
	    
}
