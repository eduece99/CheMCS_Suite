package org.cisrg.hyperstructures;

import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AnyOrderQueryBond;

public class AnyOrderQueryBondNotNull extends AnyOrderQueryBond {

	 private static final long serialVersionUID = -827101570208878645L;
	    
	    public AnyOrderQueryBondNotNull(IChemObjectBuilder builder) {
	    	super(builder);
	    }
	    
	    /**
	     * Creates a new instance
	     *
	     * @param atom1
	     * @param atom2
	     */
	    public AnyOrderQueryBondNotNull(IQueryAtom atom1, IQueryAtom atom2, IBond.Order order, IChemObjectBuilder builder) {
	        super(atom1, atom2, order, builder);
	    }
	    
	    
	    public boolean matches(IBond bond) {
	            return bond != null; // any bond order is fine as long as it is not null
	        }
	    
	    public String toString() {
	    	return "AnyOrderQueryBondNotNull()";
	    }
	
}
