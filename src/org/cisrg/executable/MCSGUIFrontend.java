package org.cisrg.executable;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.DefaultComboBoxModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SpinnerModel;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.border.Border;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JFileChooser;


import org.cisrg.knime.GraphSimilarityNodeModel.AggregationMethod;
import org.cisrg.knime.GraphSimilarityNodeModel.ReturnType;
import org.cisrg.knime.GraphSimilarityNodeModel.SimilarityType;
import org.cisrg.mapping.ConvenienceTools;
import org.cisrg.mapping.ExtendedAlgorithm;
import org.cisrg.mapping.SimilarityComparator;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import picocli.CommandLine.Option;

import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;



public class MCSGUIFrontend extends Frame {
	
	
	private void addMoleculesToPanel() {
		
		
		algorithmName = mappingAlgorithmComponent.getSelectedItem().toString();
	    topologicalDistanceLimit = Integer.parseInt( topoDistanceLimitComponent.getText() );
	    ringHeuristics = ringHeuristicsComponent.isSelected();
	    raymondHeuristics = raymondHeuristicsComponent.isSelected() ;
	    timeLimit = Integer.parseInt( timeLimitComponent.getText() );
		
		ExtendedAlgorithm algorithm = ExtendedAlgorithm.valueOf(algorithmName);
		
		simHub = new SimilarityComparator(
				null, bondWeightFlag, ghostSubstructures, algorithm, 
				raymondHeuristics, ringHeuristics, topologicalDistanceLimit, timeLimit, false
		);
		
		
		
		Image bi = null;
		DepictionGenerator dptgen = new DepictionGenerator();
		dptgen = dptgen.withSize(molWidth, molHeight)             
			      .withMolTitle()
			      .withTitleColor(Color.DARK_GRAY); 
		
		GridBagConstraints gbc = new GridBagConstraints();
	    gbc.gridx = 0;
        gbc.gridy = 0;
        
        ImageIcon[] tableColumns = new ImageIcon[refMols.size()] ;
        Object[][] tableData = new Object[dbMols.size()][refMols.size()] ;
		
        // have ref molecules as columns, and db molecules as rows
		for( int r = 0; r < refMols.size(); r++ ) {
			
			try {
				Image refImg = dptgen.depict( refMols.get(r) ).toImg();
				tableColumns[r] = new ImageIcon( refImg );
				//molPanel.add( new JLabel( new ImageIcon(refImg) ), gbc );
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			for( int d = 0; d < dbMols.size(); d++ ) {
			
				simHub.calculateSimilarity(refMols.get(r), dbMols.get(d));
				
				//IAtomContainer mol = mols.get(mi);
				
				try {
					  /*
					  if( mi == 0 ) {
						  dptgen = dptgen.withHighlight( simHub.bondMaps.get(0).keySet() , Color.RED).withOuterGlowHighlight(2.0) ;
					  } else {
						  dptgen = dptgen.withHighlight( simHub.bondMaps.get(0).values() , Color.RED).withOuterGlowHighlight(2.0) ;
					  }
					  */ 
					 
					  dptgen = dptgen.withTitleScale(1.2);
					  dptgen = dptgen.withHighlight( simHub.bondMaps.get(0).keySet() , mcsColours[r] ).withOuterGlowHighlight(2.0) ;
			    	  dptgen = dptgen.withHighlight( simHub.bondMaps.get(0).values() , mcsColours[r] ).withOuterGlowHighlight(2.0) ;
					  List<IAtomContainer> molPair = new ArrayList<IAtomContainer>(2);
					  molPair.add(refMols.get(r) );
					  molPair.add(dbMols.get(d) );
					  
					  String caption = "";
					  if( similarityTypeComponent.getSelectedItem() == SimilarityType.Tanimoto ) {
						  caption =  String.valueOf(simHub.tanimoto);
					  } else if( similarityTypeComponent.getSelectedItem() == SimilarityType.Tversky ) {
						  caption =  String.valueOf(simHub.tversky);
					  } else if( similarityTypeComponent.getSelectedItem() == SimilarityType.MCSSize ) {
						  caption =  String.valueOf(simHub.bondMaps.get(0).size() );
					  } else if( similarityTypeComponent.getSelectedItem() == SimilarityType.MCSTime ) {
						  caption =  String.valueOf(simHub.mcsExecTime);
					  } else if( similarityTypeComponent.getSelectedItem() == SimilarityType.FragmentSizes ) {
						  caption =  String.valueOf( Arrays.stream(simHub.fragmentSizes).boxed().toList() );
					  }
					  refMols.get(r).setProperty(CDKConstants.TITLE, caption );
					  
			    	  bi = dptgen.depict( molPair ).toImg();
					  
					  //molPanel.add( new JLabel( new ImageIcon(bi) ), gbc );
					  
					  tableData[d][r] = new ImageIcon( bi );
					  gbc.gridy++;
								
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			  
		        
			}
			
			gbc.gridx++;
		}
		
		molPanel.getViewport().add( createMolTable(tableData, tableColumns));
		//molPanel.repaint();
		//molPanel.validate();
	}
	
	
	class FileButtonListener implements ActionListener {

		JButton buttonComponent;

		public FileButtonListener( JButton bc ) {
			super();
			
			this.buttonComponent = bc;
		}
		
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			
			currentButton = buttonComponent;
			fc.showOpenDialog(buttonComponent);
			
		}
		
	}
	
	class LoadMoleculesButtonListener implements ActionListener {

		public LoadMoleculesButtonListener() {
			super();
		}
		
		@Override
		public void actionPerformed(ActionEvent ev) {
			// TODO Auto-generated method stub
			
			// get molecules
	        ArrayList<IAtomContainer> compounds = null;
			try {
				int refLim = (Integer) refMolLimComponent.getValue();
				int dbLim = (Integer) dbMolLimComponent.getValue();
				
				refMols = ConvenienceTools.getQueryMolecules( new File( refFileComponent.getText() ), null, true, refLim ); 
				dbMols = ConvenienceTools.getQueryMolecules( new File( dbFileComponent.getText() ), null, true, dbLim );
				//System.out.println( "compound file size - " + compounds.size() );
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				
				JOptionPane.showMessageDialog(null,"Input molecule files do not exist.  Check your file paths","Error",1);
				
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				
				JOptionPane.showMessageDialog(null,"Cannot read input molecules.  Only .smi and .sdf Files are supported","Error",1);
				
				e.printStackTrace();
			}
			
			
			// remove all molecules
			molPanel.getViewport().removeAll();
			
			
			addMoleculesToPanel( );
			
		}
		
	}
	
	class FileChooserListener implements ActionListener {
		
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			
			if( currentButton == refFileButtonComponent ) {
				refFileComponent.setText( fc.getSelectedFile().getAbsolutePath() );
			} else if( currentButton == dbFileButtonComponent ) {
				dbFileComponent.setText( fc.getSelectedFile().getAbsolutePath() );
			} else {
				System.out.println("No file path saved." + currentButton );
				
				
			}
		}
		
	}

	
	JComboBox mappingAlgorithmComponent;
	JComboBox aggregationMethodComponent;
	JComboBox returnTypeComponent;
	JComboBox similarityTypeComponent;
	
	JTextField topoDistanceLimitComponent;
	JTextField timeLimitComponent;
	JTextField refFileComponent;
	JTextField dbFileComponent;
	JSpinner refMolLimComponent;
	JSpinner dbMolLimComponent;
	JButton refFileButtonComponent;
	JButton dbFileButtonComponent;
	JButton loadMoleculesButtonComponent;
	
	JButton currentButton;
	
	JScrollPane molPanel;
	Color[] mcsColours = {
	  Color.RED,
	  Color.BLUE,
	  Color.MAGENTA,
	  Color.CYAN,
	  Color.GREEN
	};
	
	JCheckBox outputMCS;
	JCheckBox raymondHeuristicsComponent;
	JCheckBox ringHeuristicsComponent;
	JCheckBox detailedGhostInfoComponent;
	
	List<IAtomContainer> refMols, dbMols;
	
	final JFileChooser fc;
	
	SimilarityComparator simHub = null;
	private String algorithmName = "Depolli_dMCES";
    private int topologicalDistanceLimit = 4;
    private boolean ringHeuristics = false;
    private boolean raymondHeuristics = false;
    private int timeLimit = 10000;
	private boolean ghostSubstructures = false;
	private boolean bondWeightFlag = false;
	
	private int molHeight = 175, molWidth = 300;
	
	
	public static ExtendedAlgorithm[] MappingAlgorithmNames = new ExtendedAlgorithm[] {
			ExtendedAlgorithm.DEFAULT,
			ExtendedAlgorithm.CDKMCS,
			ExtendedAlgorithm.VFLibMCS,
			ExtendedAlgorithm.MCSPlus,
			ExtendedAlgorithm.ChemAxon_cMCES,
			ExtendedAlgorithm.ChemAxon_dMCES,
			ExtendedAlgorithm.consR_dMCES,
			ExtendedAlgorithm.BK_dMCES,
			ExtendedAlgorithm.BK_cMCES,
			ExtendedAlgorithm.CP_dMCES,
			ExtendedAlgorithm.RASCAL_dMCES,
			ExtendedAlgorithm.Depolli_dMCES,
			ExtendedAlgorithm.kCombu_dMCES,
			ExtendedAlgorithm.kCombu_cMCES,
			ExtendedAlgorithm.fMCS
		};
		

	/** Enum for the different aggregation methods. */
	public enum AggregationMethod {
		Minimum, Maximum, Average, Matrix
	}
	
	
	/** Enum for the similarity measure types. */
	public enum SimilarityType { 
		Tanimoto, Tversky, MCSSize, MCSTime, FragmentSizes
	}

	public MCSGUIFrontend() {
		
		/*
		refMolColumnComponent = new ColumnSelectionPanel( 
    			(Border) null, 
    			//new DataValueColumnFilter( StringValue.class ) ,
    			new DataValueColumnFilter( CDKValue.class, SmartsValue.class ) ,
    			//new DataValueColumnFilter( accepted ) ,
    			false  // dummy column if true
    	);
        
        databaseMolColumnComponent = new ColumnSelectionPanel( 
    			(Border) null, 
    			new DataValueColumnFilter( CDKValue.class ) ,
    			//new DataValueColumnFilter( CDKNodeUtils.ACCEPTED_VALUE_CLASSES ) ,
    			false  // dummy column if true
    	);
       */
		
		
		ExtendedAlgorithm algorithm = ExtendedAlgorithm.valueOf(algorithmName);
		
		simHub = new SimilarityComparator(
				null, bondWeightFlag, ghostSubstructures, algorithm, 
				raymondHeuristics, ringHeuristics, topologicalDistanceLimit, timeLimit, false
		);
		
		refFileComponent = new JTextField("input file path", 200);  
		refFileButtonComponent = new JButton("Browse");
		refFileButtonComponent.addActionListener( new FileButtonListener(refFileButtonComponent) );
		
		SpinnerModel refMolLimModel = new SpinnerNumberModel(5, //initial value
                1, //min
                5, //max
                1);//step
		
		SpinnerModel dbMolLimModel = new SpinnerNumberModel(100, //initial value
                1, //min
                100, //max
                1);//step


		refMolLimComponent = new JSpinner(refMolLimModel);
		dbMolLimComponent = new JSpinner(dbMolLimModel);
		
		dbFileComponent = new JTextField("input file path", 200); 
		dbFileButtonComponent = new JButton("Browse");
		dbFileButtonComponent.addActionListener( new FileButtonListener(dbFileButtonComponent) );
		
		
		loadMoleculesButtonComponent = new JButton("Load Molecules");
		loadMoleculesButtonComponent.addActionListener( new LoadMoleculesButtonListener() );
		
		
		fc = new JFileChooser();
		fc.addActionListener( new FileChooserListener( ) );
		fc.addActionListener( new FileChooserListener( ) );
       
        mappingAlgorithmComponent = new JComboBox( );
        mappingAlgorithmComponent.setModel( new DefaultComboBoxModel( MappingAlgorithmNames  ) );
        mappingAlgorithmComponent.setSelectedItem(ExtendedAlgorithm.Depolli_dMCES);
        
        aggregationMethodComponent = new JComboBox();
        aggregationMethodComponent.setModel( new DefaultComboBoxModel( AggregationMethod.values() ) );
        
        similarityTypeComponent = new JComboBox();
        similarityTypeComponent.setModel( new DefaultComboBoxModel( SimilarityType.values() ) );
        similarityTypeComponent.setSelectedIndex(0);
        
        outputMCS = new JCheckBox();
        outputMCS.setSelected(false);
        
        topoDistanceLimitComponent = new JTextField(2);  // 2 columns
        topoDistanceLimitComponent.setText("-1");
        
        timeLimitComponent = new JTextField(6);  // 6 columns
        timeLimitComponent.setText("10000");
        
        raymondHeuristicsComponent = new JCheckBox();
        raymondHeuristicsComponent.setSelected(false);
        
        ringHeuristicsComponent = new JCheckBox();
        ringHeuristicsComponent.setSelected(false);
        
        detailedGhostInfoComponent = new JCheckBox(); 
        detailedGhostInfoComponent.setSelected(false);
        		
        
        GridBagConstraints gbc = new GridBagConstraints();
        
        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new GridBagLayout() );
        mainPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        
        molPanel = new JScrollPane() ;
        molPanel.setMinimumSize( new Dimension(250, 500));
        
        //molPanel.setLayout(new GridBagLayout() );
        //molPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        
        JTabbedPane tp = new JTabbedPane();  
        //tp.setLayout( new GridBagLayout()  );
        tp.setBounds(50,50,200,200);  
        
        JPanel ioPanel = new JPanel();  
        ioPanel.setLayout(new GridBagLayout() );
        
        JPanel mcsSettingsPanel = new JPanel();  
        mcsSettingsPanel.setLayout(new GridBagLayout() );
        
        
        
        tp.add("IO", ioPanel);
        tp.add("MCS", mcsSettingsPanel);
        
        
        
        //addTab("Settings", mainPanel);
        
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weightx = 0.2;
        gbc.gridx = 0;
        gbc.gridy = 0;
        ioPanel.add(new JLabel("reference molecules: "), gbc);
        gbc.gridy++;
        
        ioPanel.add(new JLabel("database molecules: "), gbc);
        gbc.gridy++;
        
        gbc.gridx = 1;
        gbc.gridy = 0;
        
        gbc.fill = GridBagConstraints.BOTH ;
        gbc.weightx = 0.8;
        ioPanel.add( refFileComponent, gbc);
        gbc.gridy++;
        
        ioPanel.add( dbFileComponent, gbc);
        gbc.gridy++;
        
        gbc.weightx = 0.02;
        gbc.gridx = 2;
        gbc.gridy = 0;
        ioPanel.add( refMolLimComponent, gbc);
        gbc.gridy++;
        
        ioPanel.add( dbMolLimComponent, gbc);
        gbc.gridy++;
        
        
        gbc.weightx = 0.05;
        gbc.gridx = 3;
        gbc.gridy = 0;
        ioPanel.add( refFileButtonComponent, gbc);
        gbc.gridy++;
        
        ioPanel.add( dbFileButtonComponent, gbc);
        gbc.gridy++;
        
        
        gbc.gridx = 0;
        ioPanel.add( loadMoleculesButtonComponent, gbc);
        gbc.gridy++;
        
        
        gbc.gridy = 0;
        mcsSettingsPanel.add(new JLabel("search algorithm: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("aggregation method: "), gbc);
        gbc.gridy++;
        
        //mainPanel.add(new JLabel("return method: "), gbc);
        //gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("similarity equation: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("Modular Product heuristics: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("Modular Product Ring heuristics: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("topological distance limit: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("expansion time limit: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("output MCS: "), gbc);
        gbc.gridy++;
        
        //mcsSettingsPanel.add(new JLabel("output ghost SMARTS: "), gbc);
        //gbc.gridy++;
        
        
        
        gbc.gridx = 1;
        gbc.gridy = 0;
             
        mcsSettingsPanel.add( mappingAlgorithmComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( aggregationMethodComponent, gbc );
        gbc.gridy++;
        
        //mainPanel.add( returnTypeComponent, gbc );
        //gbc.gridy++;
        
        mcsSettingsPanel.add( similarityTypeComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( raymondHeuristicsComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( ringHeuristicsComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( topoDistanceLimitComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( timeLimitComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( outputMCS, gbc );
        gbc.gridy++;
        
        //mcsSettingsPanel.add( detailedGhostInfoComponent, gbc );
        //gbc.gridy++;
		
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent windowEvent){
               System.exit(0);
            }        
         });    
        
        
        IChemObjectBuilder bldr   = SilentChemObjectBuilder.getInstance();
        SmilesParser       smipar = new SmilesParser(bldr);
        
        refMols = new ArrayList<IAtomContainer>(1);
        dbMols = new ArrayList<IAtomContainer>(1);
        
        IAtomContainer mol = null;
        BufferedImage bi = null;
		try {
			mol = smipar.parseSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine");
			mol.setProperty(CDKConstants.TITLE, "caffeine"); // title already set from input!
			refMols.add(mol);
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			mol = smipar.parseSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)CC methylcaffeine");
			//mol.setProperty(CDKConstants.TITLE, "caffeine"); // title already set from input!
			dbMols.add(mol);
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        

        DepictionGenerator dptgen = new DepictionGenerator();
        try {
			 bi = dptgen.withSize(molWidth, molHeight)               // px (raster) or mm (vector)
			      .withMolTitle()
			      .withTitleColor(Color.DARK_GRAY) // annotations are red by default
			      .depict(mol).toImg();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        gbc.gridx = 0;
        gbc.gridy = 0;
        //molPanel.setBackground(Color.red);
        
        Object[][] tData = new Object[][]{{ new ImageIcon(bi) }};
        ImageIcon[] tCols = new ImageIcon[]{ new ImageIcon(bi) };
        
        addMoleculesToPanel( );
        //molPanel.getViewport().add( createMolTable(tData, tCols) );
        //molPanel.add( exampleTable, gbc );
        gbc.gridy++;
              
        
        setSize(600,600);  
        setLayout( new GridBagLayout() );
        
        
        // add panels to Frame
        
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weighty = 0.8;
        gbc.fill = gbc.BOTH;
        add( molPanel , gbc );
        
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weighty = 0.01;
        gbc.fill = gbc.HORIZONTAL;
        add( mainPanel , gbc );
        
        
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.weighty = 0.09;
        add( tp , gbc );
        
        setTitle("CheMCS GUI");
		setVisible(true);  
	}
	
	
	private JTable createMolTable( Object[][] tD, ImageIcon[] tC ) {
		DefaultTableModel model = new DefaultTableModel(tD, tC)
        {
            //  Returning the Class of each column will allow different
            //  renderers to be used based on Class
            public Class getColumnClass(int column)
            {
                return getValueAt(0, column).getClass();
            }
        };
        
        JTable molTable = new JTable( model );
        //molTable.getColumnModel().getColumn(0).setCellRenderer(new DefaultTableCellRenderer());
        //exampleTable.setPreferredScrollableViewportSize(exampleTable.getPreferredSize());
        molTable.setFillsViewportHeight(true);
        molTable.setRowHeight(molHeight);
        
        JTableHeader header = molTable.getTableHeader();
        //header.setBackground(Color.yellow);
        header.setDefaultRenderer( new IconRenderer(tC) );
        
        return(molTable);
	}
	
	class IconRenderer extends DefaultTableCellRenderer {
		
		public IconRenderer( ImageIcon[] tC ) {
			super();
			
			colIcons = tC;
		}
		
		  public Component getTableCellRendererComponent(JTable table, 
		Object obj,boolean isSelected, boolean hasFocus, int row, 
		int column) {
		   
		  setBorder(UIManager.getBorder("TableHeader.cellBorder"));
		  setHorizontalAlignment(JLabel.CENTER);
		  //System.out.println( "test " + obj.getClass() +  "  " + obj  );
		  
		  //setIcon( new ImageIcon("/home/edmund/git/Java_MCS_algorithms/data/output/chembl751606_aid466_decoys_r0_d0.png"));
		  setIcon( colIcons[column] );
		  
		  return this;
		  }
		  
		  ImageIcon[] colIcons;
	}
	
	
	
	public static void main(String[] args) {
		
		try {
            // Set System L&F
			//UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			UIManager.setLookAndFeel("com.sun.java.swing.plaf.gtk.GTKLookAndFeel");
		} 
	    catch (UnsupportedLookAndFeelException e) {
	       e.printStackTrace();
	       
	    }
	    catch (ClassNotFoundException e) {
	    	e.printStackTrace();
	    }
	    catch (InstantiationException e) {
	    	e.printStackTrace();
	    }
	    catch (IllegalAccessException e) {
	    	e.printStackTrace();
	    }
		
		MCSGUIFrontend f = new MCSGUIFrontend();    

	}
	
	

}
