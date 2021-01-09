package main

import (
    "fmt"
    "github.com/jhpeng/ctQMC/models/link_rule"
    "github.com/jhpeng/ctQMC/models/insert_rule"
)


type Model struct {
    Bond2type []int
    Bond2hNspin []int
    Bond2weight []float64
    Bond2index [][]int
    LinkRule [][]int
    InsertRule []func([]int)bool
    Nsite int
    Nbond int
}


func IsingTFSquare(lx int, ly int, hx float64) Model {
    nsite := lx*ly
    nbond1 := nsite
    nbond2 := 2*nsite
    nbond  := nbond1+nbond2

    var m Model
    m.Bond2type = make([]int,nbond)
    m.Bond2hNspin = make([]int,nbond)
    m.Bond2weight = make([]float64,nbond)
    m.Bond2index = make([][]int,nbond)

    for i:=0;i<nbond1;i++ {
        m.Bond2type[i]   = 0
        m.Bond2hNspin[i] = 1
        m.Bond2weight[i] = hx*0.5
        m.Bond2index[i] = []int{i}
    }
    n:=nbond1
    for ix:=0;ix<lx;ix++ {
        for iy:=0;iy<ly;iy++ {
            i := ix+iy*lx
            j := (ix+1)%lx+iy*lx

            m.Bond2type[n]   = 1
            m.Bond2hNspin[n] = 2
            m.Bond2weight[n] = 0.5
            m.Bond2index[n]  = []int{i,j}
            n++
        }
    }
    for ix:=0;ix<lx;ix++ {
        for iy:=0;iy<ly;iy++ {
            i := ix+iy*lx
            j := ix+((iy+1)%ly)*lx

            m.Bond2type[n]   = 1
            m.Bond2hNspin[n] = 2
            m.Bond2weight[n] = 0.5
            m.Bond2index[n]  = []int{i,j}
            n++
        }
    }


    m.LinkRule = make([][]int,2)
    m.LinkRule[0] = link_rule.LinkRuleIsing1
    m.LinkRule[1] = link_rule.LinkRuleIsing2

    m.InsertRule = make([]func([]int)bool,2)
    m.InsertRule[0] = insert_rule.InsertRuleIsing1
    m.InsertRule[1] = insert_rule.InsertRuleIsing2

    m.Nsite = nsite
    m.Nbond = nbond

    return m
}

func main() {
    m := IsingTFSquare(4,4,1.5)

    fmt.Printf("%v\n",m)
}
