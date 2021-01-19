package dtype


type Id struct {
    Key uint64
    I int
}

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
    Bond int
    HNspin int
    State []int
}

type WorldLine struct {
    Table map[uint64]Vertex
    SequenceA []uint64
    SequenceB []uint64
    Nvertices int
    Flag bool                 //true -> SequenceA; false -> SequenceB
    State []int
    Last  []Id
    First []Id
    Nsite int
    Beta float64
}

type Estimator struct {
    Name string
    Sample []float64
    N int
    Block  []float64
    Bsize  int
    Nblock int
}
