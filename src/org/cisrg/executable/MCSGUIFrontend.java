package org.cisrg.executable;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.DefaultComboBoxModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.ProgressMonitor;
import javax.swing.ScrollPaneConstants;
import javax.swing.ScrollPaneLayout;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingWorker;
import javax.swing.SpinnerModel;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.border.Border;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JFileChooser;
import javax.swing.JDialog;

import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.LongColumn;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;
import tech.tablesaw.io.csv.CsvWriter;
import tech.tablesaw.io.Destination;

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
	
	

	
	JComboBox<ExtendedAlgorithm> mappingAlgorithmComponent;
	JComboBox<String> aggregationMethodComponent;
	JComboBox<String> returnTypeComponent;
	JComboBox<String> similarityTypeComponent;
	JComboBox<MolRenderSize> molSizeComponent;
	
	JTextField topoDistanceLimitComponent;
	JTextField timeLimitComponent;
	JTextField refFileComponent;
	JTextField dbFileComponent;
	JTextField exportTableFileComponent;
	JTextField exportRenderFileComponent;
	JSpinner refMolLimComponent;
	JSpinner dbMolLimComponent;
	
	JButton refFileButtonComponent;
	JButton dbFileButtonComponent;
	JButton loadMoleculesButtonComponent;
	JButton tableBrowseButtonComponent;
	JButton renderBrowseButtonComponent;
	JButton exportTableButtonComponent;
	JButton exportRenderButtonComponent;
	
	JButton currentButton;
	
	JDialog mcsProgressDialog ;
	JLabel mcsProgressLabel ;
	
	//ProgressMonitor mcsProgressMonitor;
	int maxProgress = 1000;
	
	
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
	
	/** Enum for the molecule output sizes. */
	public enum MolRenderSize { 
		Small, Medium, Large
	}
	
	
	
	
	
	
	private StringColumn refNameCol;
	private StringColumn dbNameCol;
	private StringColumn mcsSMARTSCol;
	private IntColumn mcsSizeCol;
	private DoubleColumn mcsTanimotoCol;
	private DoubleColumn mcsTverskyCol;
	private StringColumn fragmentSizesCol;
	private IntColumn mcsTimeCol;
	
	private Collection<String> tbRefNames, tbDbNames, tbMcsSMARTS, tbFragmentSizes;
	private Collection<Integer> tbMcsSizes, tbMcsTimes;
	private Collection<Double> tbMcsTanimoto, tbMcsTversky;
	private ProgressMonitor mcsProgressMonitor;
	
	
	class MCSTask extends SwingWorker<Void, Void> {
		
		private ImageIcon[] tableColumns ;
        private Object[][] tableData;
        
        public MCSTask( ImageIcon[] tc, Object[][] td) {
        	super();
        	
        	this.tableColumns = tc;
        	this.tableData = td;
        }
        
        public Object[][] getTableData() {
        	return tableData;
        }
        
        @Override
        public Void doInBackground() {
            
        	
        	Image bi = null;
    		DepictionGenerator dptgen = new DepictionGenerator();
    		dptgen = dptgen.withSize(molWidth, molHeight)             
    			      .withMolTitle()
    			      .withTitleColor(Color.DARK_GRAY); 
    		
    		int maxMolCount = refMols.size() * dbMols.size() ;

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
    			
    			for( int d = 0; d < dbMols.size() && ! mcsProgressMonitor.isCanceled() ; d++ ) {
    			
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
    					  
    					  String caption = "no caption";
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
    					  //refMols.get(r).setProperty(CDKConstants.TITLE, caption );
    					  
    					  // Table information
    					  tbRefNames.add( refMols.get(r).getTitle() );
    					  tbDbNames.add( dbMols.get(r).getTitle() );
    					  tbMcsSMARTS.add( simHub.mcsSMARTS );
    					  tbFragmentSizes.add( String.valueOf( Arrays.stream(simHub.fragmentSizes).boxed().toList() ) ); 
    					  tbMcsSizes.add( simHub.bondMaps.get(0).size() );
    					  tbMcsTimes.add( Integer.valueOf( (int) simHub.mcsExecTime ) );
    					  tbMcsTanimoto.add( simHub.tanimoto );
    					  tbMcsTversky.add( simHub.tversky );
    					  
    					  
    			    	  bi = dptgen.depict( molPair ).toImg();
    			    	  
    			    	  Graphics2D biG = (Graphics2D) bi.getGraphics();
    			    	  //biG.setFont(new Font("TimesRoman", Font.PLAIN, 12)); 
    			    	  int noImgWidth = biG.getFontMetrics().stringWidth(caption);
    			    	  biG.setPaint(Color.BLACK);
    			    	  biG.drawString(caption, (molWidth - noImgWidth) / 2, molHeight - 5);
    					  
    					  //molPanel.add( new JLabel( new ImageIcon(bi) ), gbc );
    					  
    					  tableData[d][r] = new ImageIcon( bi );
    					  
    					  
    					  float progressIndicator = ( (float) tbMcsSMARTS.size() / maxMolCount ) * maxProgress ;
    					  mcsProgressMonitor.setProgress( Math.round(progressIndicator) );
    					  
    					  mcsProgressLabel.setText("this is molecule " + (r * d) );
    					  mcsProgressDialog.validate();
    					  
    					  //System.out.println( mcsProgressMonitor.isCanceled() );
    					  
    				} catch (CDKException e) {
    					// TODO Auto-generated catch block
    					e.printStackTrace();
    				}
    			  
    		        
    			}
    			
    			//gbc.gridx++;
    		}
    		
    		
            return null;
        }

        @Override
        public void done() {
            //Toolkit.getDefaultToolkit().beep();
            //startButton.setEnabled(true);
            mcsProgressMonitor.setProgress(0);
            mcsProgressMonitor.close();
            
            System.out.println("MCS calc done " + tableData[0][0] );
    		
    		molPanel.getViewport().add( createMolTable(tableData, tableColumns));
    		
    		
    		molPanel.repaint();
    		molPanel.validate();
        }
        
    }

	
	/**
	 * Output Table
	 * 
	 * Table :
	 * - each row contains:
	 * -- ref mol name
	 * -- db mol name
	 * -- MCS SMARTS
	 * -- MCS size
	 * -- Tanimoto similarity
	 * -- Tversky similarity
	 * -- Fragment sizes in MCS
	 * -- MCS time
	 * 
	 * If you want ref mols and db mols as columns and rows respectively, use 
	 * a competent data science tool (like python pandas) to pivot accordingly.
	 */
	private void exportTableDelim() {
		refNameCol = StringColumn.create( "reference", tbRefNames );
		dbNameCol = StringColumn.create( "database", tbDbNames );
		mcsSMARTSCol = StringColumn.create( "MCS_SMARTS", tbMcsSMARTS );
		mcsSizeCol = IntColumn.create( "MCS_bondcount", tbMcsSizes.toArray(new Integer[1]) );
		mcsTanimotoCol = DoubleColumn.create( "tanimoto", tbMcsTanimoto );
		mcsTverskyCol= DoubleColumn.create( "tversky", tbMcsTversky );
		fragmentSizesCol = StringColumn.create( "fragment_sizes", tbFragmentSizes );
		mcsTimeCol = IntColumn.create( "MCS_Time", tbMcsTimes.toArray(new Integer[1]) );
		

		Table mcsPairStats =
		    Table.create("MCS Pair Stats")
		        .addColumns( refNameCol, dbNameCol, mcsSMARTSCol, mcsSizeCol, 
		        		mcsTanimotoCol, mcsTverskyCol, fragmentSizesCol, mcsTimeCol
		        );
		
		System.out.println("Table Output:");
		//System.out.println( mcsPairStats );
		CsvWriter csvWriter = new CsvWriter();
		csvWriter.write(mcsPairStats, new Destination( new File( exportTableFileComponent.getText() ) ) );
				
	}
	
	
	
	private void addMoleculesToPanel() {
		
		
		tbRefNames = new ArrayList<String>();
		tbDbNames = new ArrayList<String>();
		tbMcsSMARTS = new ArrayList<String>();
		tbFragmentSizes = new ArrayList<String>(); 
		tbMcsSizes = new ArrayList<Integer>();
		tbMcsTimes = new ArrayList<Integer>();
		tbMcsTanimoto = new ArrayList<Double>();
		tbMcsTversky = new ArrayList<Double>();
		
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
		
		
		// define molecule dimensions
		switch( (MolRenderSize) molSizeComponent.getSelectedItem() ) {
			case Small:
				molWidth = 160;
				molHeight = 100;
				break;
			case Medium:
				molWidth = 300;
				molHeight = 175;
				break;
			case Large:
				molWidth = 400;
				molHeight = 250;
				break;
			default:
				molWidth = 160;
				molHeight = 100;
				break;
		}
				
		
		
        ImageIcon[] tableColumns = new ImageIcon[refMols.size()] ;
        Object[][] tableData = new Object[dbMols.size()][refMols.size()] ;
        
        
        
        //mcsProgressDialog.setVisible(true);
		
        mcsProgressMonitor = new ProgressMonitor( this,
                "Calculating MCS pairs for input molecule(s)...",
                "", 0, maxProgress );
        
        mcsProgressMonitor.setMillisToDecideToPopup(250);
        
        MCSTask task = new MCSTask( tableColumns, tableData );
        task.execute();
        //tableData = task.getTableData();
        
       
        
        String noImg = "No Image";
        BufferedImage image = new BufferedImage( molWidth, molHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = image.createGraphics();
        int noImgWidth = g.getFontMetrics().stringWidth(noImg);
        g.drawString(noImg, (molWidth - noImgWidth) / 2, molHeight / 2);
        //cmp.printAll(g);
        
        // fill null values (in case of premature cancellation) for rendering
        for( int i = 0; i < tableData.length; i++ ) {
        	for( int j = 0; j < tableData[i].length; j++ ) {
        		if( tableData[i][j] == null ) {
        			tableData[i][j] = new ImageIcon(image);
        		}
        	}
        }
		
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
        molTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);  // display horizontal scroll bar
        //molTable.setHorizontalScrollEnabled(true);
        
        
        //molTable.setAutoResizeMode(JTable.AUTO_RESIZE_LAST_COLUMN);
        //molTable.setMinimumSize( new Dimension( 900, 600  ));
        
        JTableHeader header = molTable.getTableHeader();
        //header.setBackground(Color.yellow);
        
        header.setDefaultRenderer( new IconRenderer(tC) );
        
        
        //molTable.setPreferredSize( new Dimension( 800, 600) );
        
        // set row height
        molTable.setRowHeight(molHeight);
        
        // set column widths
        Iterator<TableColumn> colIterator = molTable.getColumnModel().getColumns().asIterator();
        while( colIterator.hasNext() ) {
        	colIterator.next().setMinWidth( molWidth );
        }
        
        return(molTable);
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
	
	class ExportTableButtonListener implements ActionListener {
		
		@Override
		public void actionPerformed(ActionEvent e) {
			exportTableDelim();
		}
		
	}
	
	
	class ExportRenderButtonListener implements ActionListener {
		
		@Override
		public void actionPerformed(ActionEvent e) {
			
			// get the table 
			Component c = molPanel.getViewport().getComponent(0) ;
			int currentWidth = c.getWidth();
			
			// automatically set size for rendering
			int tWidth = Math.max(
					currentWidth, 
					molWidth * Math.max(1, refMols.size() )
			);
			int tHeight = molHeight * Math.max(1, dbMols.size() );
			
			c.setSize( tWidth, tHeight );
			
			//System.out.println( c + " " + c.size() );
			
			outputComponentImage( c );
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
			} else if( currentButton == tableBrowseButtonComponent ) {
				exportTableFileComponent.setText( fc.getSelectedFile().getAbsolutePath() );
			} else if( currentButton == renderBrowseButtonComponent ) {
				exportRenderFileComponent.setText( fc.getSelectedFile().getAbsolutePath() );
			} else {
				System.out.println("No file path saved." + currentButton );
			}
			
			
		}
		
	}
	

	private void outputComponentImage( Component cmp ) {
		
		BufferedImage image = new BufferedImage(cmp.getWidth(), cmp.getHeight(), BufferedImage.TYPE_INT_RGB);
        Graphics2D g = image.createGraphics();
        
        cmp.printAll(g);
        g.dispose();
        try { 
            ImageIO.write(image, "png", new File( exportRenderFileComponent.getText() ) ); 
        } catch (IOException e) {
            e.printStackTrace();
        }
        
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
		  setSize(molWidth, molHeight);
		  setMaximumSize( new Dimension( molWidth, molHeight ) );
		  
		  return this;
		}
		  
		ImageIcon[] colIcons;
	}
	


	public MCSGUIFrontend() {
		
		
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

		mcsProgressDialog = new JDialog(this, 
	            "Click a button",
	            false);
		mcsProgressDialog.setDefaultCloseOperation(
			    JDialog.HIDE_ON_CLOSE);
		
		mcsProgressLabel = new JLabel("placeholder");
		mcsProgressDialog.add( mcsProgressLabel );
		mcsProgressDialog.pack();
		//mcsProgressDialog.setVisible(true);
		
		
		

		refMolLimComponent = new JSpinner(refMolLimModel);
		dbMolLimComponent = new JSpinner(dbMolLimModel);
		
		dbFileComponent = new JTextField("input file path", 200); 
		dbFileButtonComponent = new JButton("Browse");
		dbFileButtonComponent.addActionListener( new FileButtonListener(dbFileButtonComponent) );
		
		
		loadMoleculesButtonComponent = new JButton("Load Molecules");
		loadMoleculesButtonComponent.addActionListener( new LoadMoleculesButtonListener() );
		
		exportTableFileComponent = new JTextField("table file path", 200); 
		tableBrowseButtonComponent = new JButton("Browse");
		tableBrowseButtonComponent.addActionListener( new FileButtonListener(tableBrowseButtonComponent) );
		
		exportRenderFileComponent = new JTextField("render file path", 200); 
		renderBrowseButtonComponent = new JButton("Browse");
		renderBrowseButtonComponent.addActionListener( new FileButtonListener(renderBrowseButtonComponent) );
		
		
		fc = new JFileChooser();
		fc.addActionListener( new FileChooserListener( ) );
		fc.addActionListener( new FileChooserListener( ) );
       
        mappingAlgorithmComponent = new JComboBox<ExtendedAlgorithm>( );
        mappingAlgorithmComponent.setModel( new DefaultComboBoxModel<ExtendedAlgorithm>( MappingAlgorithmNames  ) );
        mappingAlgorithmComponent.setSelectedItem(ExtendedAlgorithm.Depolli_dMCES);
        
        aggregationMethodComponent = new JComboBox<String>();
        aggregationMethodComponent.setModel( new DefaultComboBoxModel( AggregationMethod.values() ) );
        
        similarityTypeComponent = new JComboBox<String>();
        similarityTypeComponent.setModel( new DefaultComboBoxModel( SimilarityType.values() ) );
        similarityTypeComponent.setSelectedIndex(0);
        
        molSizeComponent = new JComboBox<MolRenderSize>();
        molSizeComponent.setModel( new DefaultComboBoxModel<MolRenderSize>( MolRenderSize.values() ) );
        molSizeComponent.setSelectedIndex(1);
        
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
        
        exportTableButtonComponent = new JButton("Export as CSV");
        exportTableButtonComponent.addActionListener( new ExportTableButtonListener() );		
        
        exportRenderButtonComponent = new JButton("Render to...");
        exportRenderButtonComponent.addActionListener( new ExportRenderButtonListener() );		
        
        
        GridBagConstraints gbc = new GridBagConstraints();
        
        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new GridBagLayout() );
        mainPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        
        molPanel = new JScrollPane();
        molPanel.setMinimumSize( new Dimension(750, 500));

        //molPanel.setLayout(new GridBagLayout() );
        //molPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        
        JTabbedPane tp = new JTabbedPane();  
        //tp.setLayout( new GridBagLayout()  );
        //tp.setBounds(50,50,200,200);  
        
        JPanel inputPanel = new JPanel();  
        inputPanel.setLayout(new GridBagLayout() );
        
        JPanel outputPanel = new JPanel();
        outputPanel.setLayout(new GridBagLayout() );
        
        JPanel mcsSettingsPanel = new JPanel();  
        mcsSettingsPanel.setLayout(new GridBagLayout() );
        
        
        
        tp.add("Input", inputPanel);
        tp.add("Output", outputPanel);
        tp.add("MCS", mcsSettingsPanel);
        
        
        
        //addTab("Settings", mainPanel);
        
        // Input panel settings
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weightx = 0.2;
        gbc.gridx = 0;
        gbc.gridy = 0;
        inputPanel.add(new JLabel("reference molecules: "), gbc);
        gbc.gridy++;
        
        inputPanel.add(new JLabel("database molecules: "), gbc);
        gbc.gridy++;
        
        gbc.gridx = 1;
        gbc.gridy = 0;
        
        gbc.fill = GridBagConstraints.BOTH ;
        gbc.weightx = 0.8;
        inputPanel.add( refFileComponent, gbc);
        gbc.gridy++;
        
        inputPanel.add( dbFileComponent, gbc);
        gbc.gridy++;
        
        gbc.weightx = 0.02;
        gbc.gridx = 2;
        gbc.gridy = 0;
        inputPanel.add( refMolLimComponent, gbc);
        gbc.gridy++;
        
        inputPanel.add( dbMolLimComponent, gbc);
        gbc.gridy++;
        
        
        gbc.weightx = 0.05;
        gbc.gridx = 3;
        gbc.gridy = 0;
        inputPanel.add( refFileButtonComponent, gbc);
        gbc.gridy++;
        
        inputPanel.add( dbFileButtonComponent, gbc);
        gbc.gridy++;
        
        
        gbc.gridx = 0;
        inputPanel.add( loadMoleculesButtonComponent, gbc);
        gbc.gridy++;
        
        
        
        
        // Output panel settings
        gbc.gridx = 0;
        gbc.gridy = 0;
        outputPanel.add(new JLabel("Molecule Size: "), gbc);
        gbc.gridy++;
        
        outputPanel.add(new JLabel("aggregation method: "), gbc);
        gbc.gridy++;
        
        outputPanel.add(new JLabel("similarity equation: "), gbc);
        gbc.gridy++;
        
        gbc.weightx = 0.2;
        outputPanel.add(exportRenderButtonComponent, gbc);
        gbc.gridy++;
        
        outputPanel.add(exportTableButtonComponent, gbc);
        gbc.gridy++;
        
        
        gbc.gridx = 1;
        gbc.gridy = 0;
        
        outputPanel.add( molSizeComponent, gbc);
        gbc.gridy++;
        
        outputPanel.add( aggregationMethodComponent, gbc );
        gbc.gridy++;
        
        outputPanel.add( similarityTypeComponent, gbc );
        gbc.gridy++;
        
        gbc.weightx = 0.6;
        outputPanel.add( exportRenderFileComponent, gbc );
        gbc.gridy++;       
        
        outputPanel.add( exportTableFileComponent, gbc );
        gbc.gridy--;
        
        gbc.gridx = 2;
        gbc.weightx = 0.2;
        outputPanel.add( renderBrowseButtonComponent, gbc );
        gbc.gridy++;
        outputPanel.add( tableBrowseButtonComponent, gbc );
        
        
        // MCS settings panel
        gbc.gridx = 0;
        gbc.gridy = 0;
        mcsSettingsPanel.add(new JLabel("search algorithm: "), gbc);
        gbc.gridy++;
        
        
        mcsSettingsPanel.add(new JLabel("Modular Product heuristics: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("Modular Product Ring heuristics: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("topological distance limit: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("expansion time limit: "), gbc);
        gbc.gridy++;
        
        
        
        //mcsSettingsPanel.add(new JLabel("output ghost SMARTS: "), gbc);
        //gbc.gridy++;
        
        
        
        gbc.gridx = 1;
        gbc.gridy = 0;
             
        mcsSettingsPanel.add( mappingAlgorithmComponent, gbc );
        gbc.gridy++;
        
        
        //mainPanel.add( returnTypeComponent, gbc );
        //gbc.gridy++;

        
        mcsSettingsPanel.add( raymondHeuristicsComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( ringHeuristicsComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( topoDistanceLimitComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( timeLimitComponent, gbc );
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
              
        
        setSize(600,800);  
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
		
		//outputComponentImage( molPanel );
		
		
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
