package org.cisrg.ambit;

import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.matchers.smarts.SMARTSBond;

public class DoubleBondAromaticityNotSpecified extends SMARTSBond
{
	private static final long serialVersionUID = -9341431640874352L;
	
	public DoubleBondAromaticityNotSpecified( IChemObjectBuilder builder ) {
    	super( builder );
    }
	
	public boolean matches(IBond bond) 
	{	
		if (bond.getOrder() == IBond.Order.DOUBLE)
			return(true);
					
		return false;
    };

}
