package org.cisrg.ambit;

import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.matchers.smarts.SMARTSBond;

public class TripleBondAromaticityNotSpecified extends SMARTSBond
{
	private static final long serialVersionUID = -433358060406382L;
	
	public TripleBondAromaticityNotSpecified( IChemObjectBuilder builder ) {
    	super( builder );
    }
	
	public boolean matches(IBond bond) 
	{	
		if (bond.getOrder() == IBond.Order.TRIPLE)
			return(true);
					
		return false;
    };

}
