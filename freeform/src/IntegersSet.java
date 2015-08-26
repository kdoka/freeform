
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;

public class IntegersSet extends HashSet
{
	private static final long serialVersionUID = -5429260598214037195L;
	
	public IntegersSet()
	{
		super();
	}
	
	public IntegersSet(int size)
	{
		super(size);
	}
	
	public IntegersSet(Collection c)
	{
		super(c);
	}
	
	public void add(int value)
	{
		add(new Integer(value));
	}
	
	public boolean contains(int value)
	{
		return contains(new Integer(value));
	}
	
	//////////////
	public int getMax(){
		return Collections.max(this);
	}
	public int getMin(){
		return Collections.min(this);
	}
}