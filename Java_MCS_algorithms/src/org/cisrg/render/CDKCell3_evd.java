/*
 * Copyright (C) 2003 - 2013 University of Konstanz, Germany and KNIME GmbH, Konstanz, Germany Website:
 * http://www.knime.org; Email: contact@knime.org
 * 
 * This file is part of the KNIME CDK plugin.
 * 
 * The KNIME CDK plugin is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
 * General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 * 
 * The KNIME CDK plugin is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along with the plugin. If not, see
 * <http://www.gnu.org/licenses/>.
 */
package org.cisrg.render;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringWriter;
import java.util.zip.GZIPInputStream;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import org.knime.chem.types.SdfValue;
import org.knime.chem.types.SmilesValue;
import org.knime.core.data.AdapterValue;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;
import org.knime.core.data.DataValue;
import org.knime.core.data.StringValue;
import org.knime.core.node.NodeLogger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.cml.CmlKnimeCore;
import org.openscience.cdk.knime.commons.CDKNodeUtils;
import org.openscience.cdk.knime.type.CDKAdapterCell;
import org.openscience.cdk.knime.type.CDKValue;
import org.openscience.cdk.layout.LayoutHelper;
import org.openscience.cdk.silent.ChemFile; 
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.stereo.StereoElementFactory;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

/**
 * Smiles {@link DataCell} holding a string as internal representation.
 * 
 * @author Bernd Wiswedel, University of Konstanz
 * @author Stephan Beisken, European Bioinformatics Institute
 * 
 * 
 * XXX  NOTE by Edmund Duesbury 20/3/2015 - I have modified this to accept colours from bonds via CDKConstants.ANNOTATIONS rather than the 
 * original method of only colouring the atoms, then colouring the bonds in between. 
 */
public final class CDKCell3_evd extends DataCell implements CDKValue, SmilesValue, SdfValue, StringValue {

	/**
	 * Convenience access member for <code>DataType.getType(CDKCell)</code>.
	 * 
	 * @see DataType#getType(Class)
	 */
	public static final DataType TYPE = DataType.getType(CDKCell3_evd.class);

	/**
	 * Returns the preferred value class of this cell implementation. This
	 * method is called per reflection to determine which is the preferred
	 * renderer, comparator, etc.
	 * 
	 * @return CDKValue.class
	 */
	public static final Class<? extends DataValue> getPreferredValueClass() {
		return CDKValue.class;
	}

	private static final NodeLogger LOGGER = NodeLogger.getLogger(CDKCell3_evd.class);

	/**
	 * Name of the data column spec property that indicates if the molecules in
	 * the column have 2D coordinates.
	 */
	public static final String COORD2D_AVAILABLE = "2D coordinates available";

	/**
	 * Name of the data column spec property that indicates if the molecules in
	 * the column have 3D coordinates.
	 */
	public static final String COORD3D_AVAILABLE = "3D coordinates available";

	/**
	 * Static instance of the serializer.
	 */
	private static final DataCellSerializer<CDKCell3_evd> SERIALIZER = new CDKSerializer();

	/**
	 * Returns the factory to read/write DataCells of this class from/to a
	 * DataInput/DataOutput. This method is called via reflection.
	 * 
	 * @return A serializer for reading/writing cells of this kind.
	 * @see DataCell
	 */
	public static final DataCellSerializer<CDKCell3_evd> getCellSerializer() {
		return SERIALIZER;
	}

	/**
	 * The visual representation for this CDK cell.
	 */
	private final String smiles;
	/**
	 * The hash code.
	 */
	private final long hash;
	/**
	 * The aux vector: 2d|3d;coords;atomIndex-color
	 */
	private final byte[] auxBytes;

	/**
	 * Creates a new DataCell containing the atom container.
	 * 
	 * @param atomContainer an atom container
	 * @return a new data cell containing the atom container
	 */
	public static DataCell createCDKCell(final IAtomContainer atomContainer) {
		return new CDKAdapterCell(new CDKCell3_evd(atomContainer));
	}

	/**
	 * Creates a new DataCell containing the atom container.
	 * 
	 * @param atomContainer an atom container
	 * @return a new data cell containing the atom container
	 */
	public static DataCell createCDKCell(final DataCell source, final IAtomContainer atomContainer) {
		return new CDKAdapterCell((AdapterValue) source, new CDKCell3_evd(atomContainer));
	}

	/**
	 * Creates new CDK cell.
	 * 
	 * @param atomContainer the CDK atom container
	 */
	public CDKCell3_evd(final IAtomContainer atomContainer) {

		int[] seq = new int[atomContainer.getAtomCount()];

		smiles = CDKNodeUtils.calculateSmiles(atomContainer, seq);

		if (smiles.length() == 0) { // should never happen
			hash = -1;
			auxBytes = new byte[0];
		} else {
			int[] aux = new int[seq.length];
			for (int v = 0; v < seq.length; v++) {
				aux[seq[v]] = v;
			}

			hash = CDKNodeUtils.calculateSimpleHash(atomContainer);
			auxBytes = toByte(atomContainer, aux);
		}
	}

	
	protected final int multiplier_2d3d = 5;
	protected final int multiplier_3d = 3;
	protected final int multiplier_2d = 2;
	protected final int numBytesPerUnit = 8;
	protected final int atomBondSeparationVal = -127;
	
	
	/**
	 * Some sort of representation where the IAtomContainer object is serialised into bytes.
	 * 
	 * The bytes contain atom coordinates, atom colour and bond colour information
	 * 
	 * @param atomContainer
	 * @param seq
	 * @return
	 */
	private byte[] toByte(final IAtomContainer atomContainer, int[] seq) {

		byte[] coords = new byte[0];
		byte[] acols = new byte[0];
		byte[] bcols = new byte[0];  // added this on for manual bond colouring

		if (GeometryTools.has2DCoordinates(atomContainer) && GeometryTools.has3DCoordinates(atomContainer)) {
			coords = new byte[atomContainer.getAtomCount() * 40 + 1];
			coords[0] = 3; // 2D & 3D
			for (int v = 0, k = 0, j = 0; v < seq.length; v++, k += 5) {
				Point2d p2 = atomContainer.getAtom(seq[v]).getPoint2d();
				System.arraycopy(toByte(p2.x), 0, coords, k * 8 + 1, 8);
				System.arraycopy(toByte(p2.y), 0, coords, k * 8 + 9, 8);
				Point3d p3 = atomContainer.getAtom(seq[v]).getPoint3d();
				System.arraycopy(toByte(p3.x), 0, coords, k * 8 + 17, 8);
				System.arraycopy(toByte(p3.y), 0, coords, k * 8 + 25, 8);
				System.arraycopy(toByte(p3.z), 0, coords, k * 8 + 33, 8);
				Integer color = atomContainer.getAtom(seq[v]).getProperty(CDKConstants.ANNOTATIONS, Integer.class);
				if (color != null) {
					byte[] tmpCols = new byte[acols.length + 8];
					System.arraycopy(acols, 0, tmpCols, 0, acols.length);
					acols = tmpCols;
					System.arraycopy(toByte(seq[v]), 0, acols, j * 8, 4);
					System.arraycopy(toByte(color), 0, acols, j++ * 8 + 4, 4);
				}
			}
		} else if (GeometryTools.has3DCoordinates(atomContainer)) {
			coords = new byte[atomContainer.getAtomCount() * 24 + 1];
			coords[0] = 1; // 3D
			for (int v = 0, k = 0, j = 0; v < seq.length; v++, k += 3) {
				Point3d p = atomContainer.getAtom(seq[v]).getPoint3d();
				System.arraycopy(toByte(p.x), 0, coords, k * 8 + 1, 8);
				System.arraycopy(toByte(p.y), 0, coords, k * 8 + 9, 8);
				System.arraycopy(toByte(p.z), 0, coords, k * 8 + 17, 8);
				Integer color = atomContainer.getAtom(seq[v]).getProperty(CDKConstants.ANNOTATIONS, Integer.class);
				if (color != null) {
					byte[] tmpCols = new byte[acols.length + 8];
					System.arraycopy(acols, 0, tmpCols, 0, acols.length);
					acols = tmpCols;
					System.arraycopy(toByte(seq[v]), 0, acols, j * 8, 4);
					System.arraycopy(toByte(color), 0, acols, j++ * 8 + 4, 4);
				}
			}
		} else if (GeometryTools.has2DCoordinates(atomContainer)) {
			coords = new byte[atomContainer.getAtomCount() * ( multiplier_2d * numBytesPerUnit ) + 1];
			coords[0] = 0; // 2D
			for (int v = 0, k = 0, j = 0; v < seq.length; v++, k += multiplier_2d) {
				Point2d p = atomContainer.getAtom(seq[v]).getPoint2d();
				System.arraycopy(toByte(p.x), 0, coords, k * numBytesPerUnit + 1, numBytesPerUnit);
				System.arraycopy(toByte(p.y), 0, coords, k * numBytesPerUnit + 9, numBytesPerUnit);
				Integer color = atomContainer.getAtom(seq[v]).getProperty(CDKConstants.ANNOTATIONS, Integer.class);
				if (color != null) {
					byte[] tmpCols = new byte[acols.length + 8];
					System.arraycopy(acols, 0, tmpCols, 0, acols.length);
					acols = tmpCols;
					System.arraycopy(toByte(seq[v]), 0, acols, j * 8, 4);
					System.arraycopy(toByte(color), 0, acols, j++ * 8 + 4, 4);
				}
			}
		}
		
		// bond colouring info, tacked on at the end
		bcols = new byte[1];
		bcols[0] = atomBondSeparationVal;  // separator from atoms
		for (int v = 0, j = 0; v < atomContainer.getBondCount(); v++) {
			Integer color = atomContainer.getBond(v).getProperty(CDKConstants.ANNOTATIONS, Integer.class);
			if (color != null) {
				byte[] tmpCols = new byte[bcols.length + numBytesPerUnit];  // appending 8 bytes to array
				System.arraycopy(bcols, 0, tmpCols, 0, bcols.length);  // copy original contents to new array
				bcols = tmpCols;
				System.arraycopy(toByte(v), 0, bcols, j * numBytesPerUnit + 1, 4);  
				System.arraycopy(toByte(color), 0, bcols, j++ * numBytesPerUnit + 5, 4);
			}
		}
		
		//System.out.println( "bcols length " + bcols.length  + " " + acols.length + " " + bcols[0] );

		byte[] aux = new byte[coords.length + acols.length + bcols.length];
		//byte[] aux = new byte[coords.length + acols.length];
		System.arraycopy(coords, 0, aux, 0, coords.length);
		System.arraycopy(acols, 0, aux, coords.length, acols.length);
		System.arraycopy(bcols, 0, aux, coords.length + acols.length, bcols.length);
		
		return aux;
	}

	
	
	public byte[] toByte(int value) {

		return new byte[] {
				(byte) ((value >> 24) & 0xff),
				(byte) ((value >> 16) & 0xff),
				(byte) ((value >> 8) & 0xff),
				(byte) ((value >> 0) & 0xff) };
	}

	private byte[] toByte(double value) {

		long raw = Double.doubleToRawLongBits(value);

		return new byte[] {
				(byte) ((raw >> 56) & 0xff),
				(byte) ((raw >> 48) & 0xff),
				(byte) ((raw >> 40) & 0xff),
				(byte) ((raw >> 32) & 0xff),
				(byte) ((raw >> 24) & 0xff),
				(byte) ((raw >> 16) & 0xff),
				(byte) ((raw >> 8) & 0xff),
				(byte) ((raw >> 0) & 0xff) };
	}

	/**
	 * Creates new CDK cell.
	 * 
	 * @param smiles the CML string
	 * @param hash the CDK hash
	 */
	public CDKCell3_evd(final String smiles, final long hash, final byte[] coordinates) {
		this.smiles = smiles;
		this.hash = hash;
		this.auxBytes = coordinates;
	}

	/**
	 * Returns the internal string value.
	 * 
	 * @return The string value.
	 */
	@Override
	public String getStringValue() {
		return getSmilesValue();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getSmilesValue() {
		return smiles;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getSdfValue() {

		IAtomContainer mol = getAtomContainerWithCoordinates();

		if (mol.getAtomCount() == 0) {
			return "";
		}

		SDFWriter sdfWriter = null;
		StringWriter stringWriter = null;

		try {
			stringWriter = new StringWriter();
			sdfWriter = new SDFWriter(stringWriter);

			if (mol != null && GeometryTools.has2DCoordinates(mol)) {
				LayoutHelper.adjustStereo(mol);
			}
			sdfWriter.write(mol);
		} catch (CDKException exception) {
			LOGGER.warn("SDfile conversion failed:", exception);
		} finally {
			try {
				sdfWriter.close();
				stringWriter.close();
			} catch (IOException exception) {
				LOGGER.warn("SDfile conversion failed:", exception);
			}
		}

		return stringWriter.toString();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public IAtomContainer getAtomContainer() {
		return getAtomContainerWithCoordinates();
	}

	/**
	 * {@inheritDoc}
	 */
	public IAtomContainer getAtomContainerWithCoordinates() {

		IAtomContainer molecule = CDKNodeUtils.getFullMolecule(smiles);
		if (molecule == null) {
			return null;
		}

		int nAtoms = molecule.getAtomCount();
		int nBonds = molecule.getBondCount();
		
		if (auxBytes.length == 0) {
			return molecule;
		}

		int f = (auxBytes[0] == 0) ? 2 : ((auxBytes[0] == 1) ? 3 : 5);

		double[] coords = new double[nAtoms * f];
		for (int i = 0; i < coords.length; i++) {
			int k = (i * 8) + 1;
			coords[i] = Double.longBitsToDouble(toLong(new byte[] {
					auxBytes[k],
					auxBytes[k + 1],
					auxBytes[k + 2],
					auxBytes[k + 3],
					auxBytes[k + 4],
					auxBytes[k + 5],
					auxBytes[k + 6],
					auxBytes[k + 7] }));
		}

		boolean hasAnnotations = (nAtoms * 8 * f + 1 == auxBytes.length) ? false : true;

		int endAtomByte = (nAtoms * 8 * f) + 1 + (nAtoms * 8);
		//System.out.println( nAtoms + " " + ((nAtoms * 8 * f) + 1) + " " + endAtomByte + " " + auxBytes.length + " " + hasAnnotations );
		int[] cols = new int[nAtoms];
		// start from first byte of the array (where atom colours begin, after coordinates)
		int i = nAtoms * 8 * f + 1;
		for ( ; i < auxBytes.length && auxBytes[i] != atomBondSeparationVal; i += 8) {
			//System.out.print( i + "," + auxBytes[i] + " " );
			int pos = toInt(new byte[] { auxBytes[i], auxBytes[i + 1], auxBytes[i + 2], auxBytes[i + 3] });
			int col = toInt(new byte[] { auxBytes[i + 4], auxBytes[i + 5], auxBytes[i + 6], auxBytes[i + 7] });
			cols[pos] = col;
		}
		
		
		// now for colouring bonds
		++i;  // increment to go past the "separator" thing between atoms and bonds
		for ( ; i < auxBytes.length; i += 8) {
			int pos = toInt(new byte[] { auxBytes[i], auxBytes[i + 1], auxBytes[i + 2], auxBytes[i + 3] });
			int col = toInt(new byte[] { auxBytes[i + 4], auxBytes[i + 5], auxBytes[i + 6], auxBytes[i + 7] });
			molecule.getBond(pos).setProperty(CDKConstants.ANNOTATIONS, col);
		}

		if (f == 2) {
			for (int v = 0, k = 0; v < nAtoms; v++, k += 2) {
				molecule.getAtom(v).setID("" + v);
				molecule.getAtom(v).setPoint2d(new Point2d(coords[k], coords[k + 1]));
				if (cols[v] != 0) {
					molecule.getAtom(v).setProperty(CDKConstants.ANNOTATIONS, cols[v]);
				}
			}
		} else if (f == 3) {
			for (int v = 0, k = 0; v < nAtoms; v++, k += 3) {
				molecule.getAtom(v).setID("" + v);
				molecule.getAtom(v).setPoint3d(new Point3d(coords[k], coords[k + 1], coords[k + 2]));
				if (cols[v] != 0) {
					molecule.getAtom(v).setProperty(CDKConstants.ANNOTATIONS, cols[v]);
				}
			}
		} else {
			for (int v = 0, k = 0; v < nAtoms; v++, k += 5) {
				molecule.getAtom(v).setID("" + v);
				molecule.getAtom(v).setPoint2d(new Point2d(coords[k], coords[k + 1]));
				molecule.getAtom(v).setPoint3d(new Point3d(coords[k + 2], coords[k + 3], coords[k + 4]));
				if (cols[v] != 0) {
					molecule.getAtom(v).setProperty(CDKConstants.ANNOTATIONS, cols[v]);
				}
			}
		}

		// automated bond colouring
		if (hasAnnotations) {
			int[] visited = new int[nAtoms];
			if (ConnectivityChecker.isConnected(molecule)) {
				IAtom atom = molecule.getAtom(0);
				visited[0] = atom.getProperty(CDKConstants.ANNOTATIONS) == null ? 1 : 2;
				//colorDfs(molecule, atom, cols, visited);
			} else {
				for (IAtomContainer mol : ConnectivityChecker.partitionIntoMolecules(molecule).atomContainers()) {
					IAtom atom = mol.getAtom(0);
					visited[Integer.parseInt(atom.getID())] = atom.getProperty(CDKConstants.ANNOTATIONS) == null ? 1
							: 2;
					//colorDfs(mol, atom, cols, visited);
				}
			}
		}

		return molecule;
	}

	private void colorDfs(IAtomContainer molecule, IAtom atom, int[] cols, int[] visited) {

		int v = Integer.parseInt(atom.getID());

		for (IAtom neighbor : molecule.getConnectedAtomsList(atom)) {

			int w = Integer.parseInt(neighbor.getID());

			if (visited[w] == 0) {
				visited[w] = neighbor.getProperty(CDKConstants.ANNOTATIONS) == null ? 1 : 2;
				if (visited[v] == 2 && visited[w] == 2) {
					molecule.getBond(atom, neighbor).setProperty(CDKConstants.ANNOTATIONS, cols[v]);
				}
				colorDfs(molecule, neighbor, cols, visited);
			} else if (visited[w] == 2) {
				if (visited[v] == 2) {
					molecule.getBond(atom, neighbor).setProperty(CDKConstants.ANNOTATIONS, cols[v]);
				}
			}
		}
	}

	public long toLong(byte[] data) {

		return (long) ((long) (0xff & data[0]) << 56 | (long) (0xff & data[1]) << 48 | (long) (0xff & data[2]) << 40
				| (long) (0xff & data[3]) << 32 | (long) (0xff & data[4]) << 24 | (long) (0xff & data[5]) << 16
				| (long) (0xff & data[6]) << 8 | (long) (0xff & data[7]) << 0);
	}

	public int toInt(byte[] data) {

		return (int) ((0xff & data[0]) << 24 | (0xff & data[1]) << 16 | (0xff & data[2]) << 8 | (0xff & data[3]) << 0);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected boolean equalsDataCell(final DataCell dc) {

		if (this == dc) {
			return true;
		}

		int hashCodeDc = ((CDKCell3_evd) dc).hashCode();
		if (dc instanceof CDKValue && this.hashCode() == hashCodeDc) {
			Long fullHashCode = CDKNodeUtils.calculateFullHash(this.getAtomContainer());
			Long fullHashCodeDc = CDKNodeUtils.calculateFullHash(((CDKCell3_evd) dc).getAtomContainer());
			return fullHashCode.hashCode() == fullHashCodeDc.hashCode();
		}

		return false;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int hashCode() {
		return ((Long) hash).hashCode();
	}

	/**
	 * Molecule hash is 64 bit
	 */
	public long hashCode64() {
		return hash;
	}

	/**
	 * Coordinates byte array
	 */
	public byte[] auxBytes() {
		return auxBytes;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		return getStringValue();
	}

	/**
	 * Factory for (de-)serializing a CDKCell.
	 */
	private static class CDKSerializer implements DataCellSerializer<CDKCell3_evd> {

		/**
		 * {@inheritDoc}
		 */
		@Override
		public void serialize(final CDKCell3_evd cell, final DataCellDataOutput out) throws IOException {
			out.writeUTF(cell.getSmilesValue());
			out.writeLong(cell.hashCode64());
			out.writeInt(cell.auxBytes().length);
			out.write(cell.auxBytes());
		}

		/**
		 * {@inheritDoc}
		 */
		@Override
		public CDKCell3_evd deserialize(final DataCellDataInput input) throws IOException {

			String blob = input.readUTF(); // either SMILES or CML
			byte[] bytes = blob.getBytes("ISO-8859-1");

			if (((bytes)[0] == (byte) (GZIPInputStream.GZIP_MAGIC))
					&& (bytes[1] == (byte) (GZIPInputStream.GZIP_MAGIC >> 8))) { // legacy CML cell
				String cml = blob;
				IAtomContainer mol = readCML(cml); // reads and uncompresses CML

				return new CDKCell3_evd(mol); // create new CDK cell
			} else { // current SMILES cell
				String smiles = blob;
				long hash64 = input.readLong();
				byte[] coords = new byte[input.readInt()];
				input.readFully(coords);

				return new CDKCell3_evd(smiles, hash64, coords);
			}
		}
	}

	/**
	 * Reads legacy compressed CML from KNIME-CDK versions 2.8 and 2.9.
	 * 
	 * @param compressedCml the compressed CML string
	 * @return an updated and fully configured CDK atom container
	 */
	private static IAtomContainer readCML(final String compressedCml) {

		IAtomContainer mol = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);

		if (compressedCml == null || compressedCml.length() == 0) {
			return mol;
		}

		try {
			GZIPInputStream gis = new GZIPInputStream(new ByteArrayInputStream(compressedCml.getBytes("ISO-8859-1")));
			CMLReader reader = new CMLReader(gis);
			reader.registerConvention(CmlKnimeCore.CONVENTION, new CmlKnimeCore());

			IChemFile chemFile = (ChemFile) reader.read(new ChemFile());
			mol = ChemFileManipulator.getAllAtomContainers(chemFile).get(0);
			mol = CDKNodeUtils.getFullMolecule(mol); // 'update' to new CDK molecules
			if (GeometryTools.has2DCoordinates(mol)) {
				StereoElementFactory stereoFactory = StereoElementFactory.using2DCoordinates(mol);
				for (IStereoElement stereoEleemnt : stereoFactory.createAll()) {
					mol.addStereoElement(stereoEleemnt);
				}
			} else if (GeometryTools.has3DCoordinates(mol)) {
				StereoElementFactory stereoFactory = StereoElementFactory.using3DCoordinates(mol);
				for (IStereoElement stereoEleemnt : stereoFactory.createAll()) {
					mol.addStereoElement(stereoEleemnt);
				}
			}

			gis = null;
			reader = null;
			chemFile = null;
		} catch (Exception exception) {
			LOGGER.warn("Deserialization of legacy CML failed.", exception);
		}

		return mol;
	}
}
