
public class HeapNode {
	int i;
	int j;
	HeapNode suc;
	public HeapNode getSuc() {
		return suc;
	}

	public void setSuc(HeapNode suc) {
		this.suc = suc;
	}

	public HeapNode getPred() {
		return pred;
	}

	public void setPred(HeapNode pred) {
		this.pred = pred;
	}

	HeapNode pred;
	double cost;

	public HeapNode(int i, int j, double cost) {
		this.i = i;
		this.j = j;
		this.cost = cost;
	}

	public int getI() {
		return i;
	}

	public int getJ() {
		return j;
	}

	public double getCost() {
		return cost;
	}
	

}
