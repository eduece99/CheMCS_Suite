package org.cisrg;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.ConcurrentModificationException;
//import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

public class BitSetExtended<E> implements Set<E> {

	
	protected BitSet elementData;
	protected int size = 0;
	protected int modCount = 0;
	
	public BitSetExtended(Collection c) {
		//elementData = new BitSet()
		
		int nBits = collectionSize( c );
		
		size = nBits;
		
		elementData = new BitSet( nBits );
		
		addAll( c );
		/*if( c instanceof BitSet ) {
			or( (BitSet) c );
		} else {
			for( Object o : c ) {
				set( o.hashCode() );
			}
		}*/

	}
	
	public BitSetExtended( BitSet bs ) {
		size = bs.length();
		elementData = bs;
	}
	
	public BitSetExtended( int s ) {
		elementData = new BitSet( s );
	}
	
	
	private int collectionSize( Collection c ) {
		
		int cSize = 0;
		
		if( c.size() == 0 )
			return 0;
		
		if( c instanceof BitSetExtended ) {
			cSize = ((BitSetExtended) c).getBitSet().size();
		} else if( c instanceof List ){
			List cList = (List) c;
			for( Object o : cList ) { cSize = Math.max( cSize, o.hashCode() ); }
			//cSize = cList.get( cList.size()-1 ).hashCode();
		} else {
			cSize = c.size();
		}
		
		return cSize;
	}
	
	
	private BitSet bitsetFromCollection( Collection c ) {
		
		if( c instanceof BitSetExtended ) {
			BitSet currentBs = 	((BitSetExtended) c).getBitSet();	
			return (BitSet) currentBs.clone();
		}
		
		BitSet other = new BitSet( collectionSize( c ) );
		
		for( Object o : c ) {
			other.set( o.hashCode() );
		}
		//System.out.println( "new bitset " + other + " " + other.length() );
		return other;
	}
	
	@Override
	public boolean add(Object e) {
		
		elementData.set( e.hashCode() );
		modCount++;
		
		return true;
	}

	// analogous to bitwise 'or' so we turn the other collection into a BitSet
	@Override
	public boolean addAll(Collection c) {
		
		int cSize = collectionSize( c );
		if( cSize == 0 )
			return false;
		
		BitSet other = bitsetFromCollection( c );
		//System.out.println( "other " + other);
		
		elementData.or( other );
		modCount++;
		
		return true;
	}

	@Override
	public boolean contains(Object o) {
		
		return elementData.get( o.hashCode() );
	}

	@Override
	public boolean containsAll(Collection c) {
		BitSet other = bitsetFromCollection( c );
		return elementData.intersects( other );
	}

	@Override
	public Iterator iterator() {
		return new Itr();
	}
	
	/*
	 * Copied this from ArrayList grepCode for Java 7-b147
	 */
	private class Itr<E> implements Iterator<E> {
        int cursor = elementData.nextSetBit(0);       // index of next element to return
        int lastRet = cursor; // index of last element returned; -1 if no such
        int expectedModCount = modCount;
        //int onBits = elementData.cardinality();
        int bsSize = elementData.length();

        public boolean hasNext() {
            return cursor >= 0 & cursor < bsSize;
        	//return (cursor != onBits || cursor >= 0);
        	//return nextSetBit(cursor+1) >= 0;
        }

        @SuppressWarnings("unchecked")
        public E next() {
            checkForComodification();
            
            int i = cursor;
            
            //if( get(i) == false )
            //System.out.println( "i " + i + " of " + onBits + " " + hasNext() );
            if( i >= bsSize || i < 0 )
            	throw new NoSuchElementException();
            
            if ( i >= bsSize )
                throw new ConcurrentModificationException();
            
            
            /*if ( i >= size() )
                throw new ConcurrentModificationException();*/
            cursor = elementData.nextSetBit( i + 1 );
            lastRet = elementData.nextSetBit(i);
            //System.out.println( "  " + cursor + " " + lastRet + " " + bsSize + " of " + hasNext() );
            return (E) new Integer(lastRet);
        }

        public void remove() {
            if (lastRet < 0)
                throw new IllegalStateException();
            checkForComodification();

            try {
                BitSetExtended.this.remove(lastRet);
                cursor = lastRet;
                lastRet = -1;
                expectedModCount = modCount;
            } catch (IndexOutOfBoundsException ex) {
                throw new ConcurrentModificationException();
            }
        }

        final void checkForComodification() {
            if (modCount != expectedModCount)
                throw new ConcurrentModificationException();
        }
    }

	@Override
	public boolean remove(Object o) {
		if( o != null ) {
			int index = o.hashCode();
			
			if( elementData.get(index) ) {
				elementData.clear( index );
				modCount++;
				return true;
			}
		}
		
		return false;
	}

	@Override
	public boolean removeAll(Collection c) {
		
		int cSize = collectionSize( c );
		if( cSize == 0 )
			return false;
		
		BitSet other = bitsetFromCollection( c );
		
		elementData.andNot( other );
		
		modCount++;
		
		return true;  // bad return type (true = Collection modified), need to fix
	}

	@Override
	public boolean retainAll(Collection c) {
		
		//int cSize = collectionSize( c );
		
		// intersection with nothing = empty set
		if( size() == 0 ) {
			elementData.clear();
			return true;
		}
		
		BitSet other = bitsetFromCollection( c );
		
		elementData.and( other );
		
		modCount++;
		
		return true;  // bad return type (true = Collection modified), need to fix
	}
	
	 

	@Override
	public Object[] toArray() {
		
		//checkInvariants();

        int numBits = elementData.cardinality() ;
        Integer[] b = new Integer[numBits];

        int i = elementData.nextSetBit(0);
        int count = 0;
        if (i != -1) {
            b[count] = i;
            for (i = elementData.nextSetBit(i+1); i >= 0; i = elementData.nextSetBit(i+1)) {
                b[++count] = i;
            }
        }

        return b;
	}

	@Override
	public Object[] toArray(Object[] a) {

		Integer[] tempAr = (Integer[]) toArray();

		if (a.length < tempAr.length) // Make a new array of a's runtime type, but my contents:
			return tempAr;

		System.arraycopy(tempAr, 0, a, 0, tempAr.length);

		if (a.length > tempAr.length)
			a[tempAr.length] = null;

		return a;
	}

	@Override
	public void clear() {
		elementData.clear();
		
		modCount++;
	}

	@Override
	public boolean isEmpty() {
		return elementData.isEmpty();
	}

	@Override
	public int size() {
		return elementData.cardinality() ;
	}
	
	
	@Override
	public String toString() {
		return elementData.toString();
	}
	
	
	public BitSet getBitSet() {
		return elementData;
	}
	
	public int getFirstBit() {
		return elementData.nextSetBit(0);
	}

	
	public static void main(String[] args) {
		
		Collection<Integer> list1 = new BitSetExtended<Integer>(10);
		Collection<Integer> list2 = new BitSetExtended<Integer>(10);
		
		list1.add(2);
		list1.add(4);
		list1.add(5);
		list1.add(7);
		
		list2.add(1);
		list2.add(5);
		list2.add(6);
		list2.add(9);
		
		/*Collection<Integer> bList1 = new BitSetExtended<Integer>(list1);
		Collection<Integer> bList2 = new BitSetExtended<Integer>(list2);*/
		Collection<Integer> bList1 = list1;
		Collection<Integer> bList2 = list2;
		
		System.out.println( "list1 " + bList1 + " " + bList1.size() );
		System.out.println( "list2 " + bList2 + " " + bList2.size() );
		//bList2.clear();
		
		bList1.addAll(bList2);
		bList1.add(12);
		System.out.println( bList1 + " " + bList1.contains(4) );
		
		for( Integer n : bList1 ) {
			System.out.print( n + " " );
		}
		System.out.println();
		
		// intersection test
		bList1.remove( ((BitSetExtended) bList1).getFirstBit() );
		System.out.println( "list1 " + bList1 + " " + bList1.size() );
		bList1.retainAll(bList2);
		System.out.println( bList1 + " " + bList2 );
		
		
		// intersection test with non bit-set
		list1.add(2);
		Collection<Integer> list3 = new ArrayList<Integer>();
		list3.add(2);
		list3.add(2);
		list3.add(4);
		list3.add(6);
		list3.add(6);
		bList1.retainAll(list3);
		System.out.println( bList1 + " " + list3 );
		
		bList1.clear();
		System.out.println( bList1 + " " + bList1.isEmpty() + " " + bList1.size() );
		
		
	}
	

}
