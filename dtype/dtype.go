package dtype

type Model struct {
    Bond2type []int
    Bond2hNspin []int
    Bond2weight []float64
    Bond2index [][]int
    LinkRule [][]int
    InsertRule []func([]int)bool
    Nsite int
    Nbond int
    Tweight float64
}

type Vertex struct {
    Key uint64
    Bond int
    HNspin int
    State []int
}

type WorldLine struct {
    SequenceA []Vertex
    SequenceB []Vertex
    Cluster   []int
    Weight    []int
    Nvertices   int
    Mnspin      int
    Flag        bool        //true -> SequenceA; false -> SequenceB
    State     []int
    Last      []int
    First     []int
    Nsite       int
    Beta        float64
}

type Estimator struct {
    Name string
    Sample []float64
    N int
    Block  []float64
    Bsize  int
    Nblock int
}
