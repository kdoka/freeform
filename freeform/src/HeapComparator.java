import java.util.Comparator;


public class HeapComparator implements Comparator {
	
	@Override
	public int compare(Object o1, Object o2) {
		if (((HeapNode) o1).getCost()<((HeapNode)o2).getCost())
	    	  return -1;
		else if (((HeapNode) o1).getCost()>((HeapNode)o2).getCost())
			return 1;
		else
	      return 0;
	}
}
