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


