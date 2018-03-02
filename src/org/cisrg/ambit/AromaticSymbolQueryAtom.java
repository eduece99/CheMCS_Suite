/*
Copyright (C) 2007-2008  

Contact: nina@acad.bg

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 2.1
of the License, or (at your option) any later version.
All we ask is that proper credit is given for our work, which includes
- but is not limited to - adding the above copyright notice to the beginning
of your source code files, and to any copyright notice that you may distribute
with programs based on this work.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA
*/
package org.cisrg.ambit;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.matchers.smarts.SMARTSAtom;

/**
 * 
 * @author Nikolay Kochev nick@uni-plovdiv.bg
 */
public class AromaticSymbolQueryAtom extends SMARTSAtom {
    
    private static final long serialVersionUID = -56452386792649224L;

    public AromaticSymbolQueryAtom( IChemObjectBuilder builder ) {
    	super( builder );
    }
    
	public boolean matches(IAtom atom) {
        if (atom.getFlag(CDKConstants.ISAROMATIC) && this.getSymbol().equals(atom.getSymbol())) 
        {	        	
        	return true;
        }
        else
        {	
        	return false;
        }	
    };

    public String toString() {
		return "AromaticSymbolQueryAtom()";
    }
}
