package org.cisrg.ambit;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.matchers.smarts.SMARTSAtom;


public class SymbolQueryAtomAromaticityNotSpecified extends SMARTSAtom {
	
	private static final long serialVersionUID = -18035318274561647L;

    public SymbolQueryAtomAromaticityNotSpecified( IChemObjectBuilder builder ) {
    	super( builder );
    }
    
	public boolean matches(IAtom atom) {
		
		if (this.getSymbol().equals(atom.getSymbol()))
        	return true;
        else
        	return false;
    };

    public String toString() {
		return "SymbolQueryAtomAromaticityNotSpecified()";
    }
}

