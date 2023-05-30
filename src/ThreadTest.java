import java.util.Random;


public class ThreadTest {
	
	public static int tCount = 0;

	public static void main( String[] args ) {
		
		
		//Thread[] threads = new Thread[3];
		final Random random = new Random();
		
		for( int t = 0; t < 18; t++ ) {
			tCount++;
			Thread tempThread = new Thread( new Runnable() {
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					
					try {
						Thread.sleep(600 + random.nextInt(2000) );
						System.out.println("Thread Complete " + tCount );
						tCount--;
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}  );
			
			tempThread.start();
			
			if( tCount > 2 ) {
				try {
					//tCount--;
					tempThread.join();
					
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				
			}
			
		}
		
	}
	
}
	
}
