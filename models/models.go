package models

import (
    "fmt"
    "github.com/jhpeng/ctQMC/models/link_rule"
    "github.com/jhpeng/ctQMC/models/insert_rule"
    . "github.com/jhpeng/ctQMC/dtype"
)


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
    m.Tweight = 0

    for i:=0;i<nbond1;i++ {
        m.Bond2type[i]   = 0
        m.Bond2hNspin[i] = 1
        m.Bond2weight[i] = hx*0.5
        m.Tweight += hx*0.5
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
            m.Tweight += 0.5
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
            m.Tweight += 0.5
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

func JQ3LadderSquare(lx int, ly int, q float64) Model {
    nsite := lx*ly
    nbond1 := 2*nsite
    nbond2 := 2*nsite
    nbond  := nbond1+nbond2

    var m Model
    m.Bond2type = make([]int,nbond)
    m.Bond2hNspin = make([]int,nbond)
    m.Bond2weight = make([]float64,nbond)
    m.Bond2index = make([][]int,nbond)
    m.Tweight = 0

    n:=0
    for iy:=0;iy<ly;iy++ {
        for ix:=0;ix<lx;ix++ {
            i := ix+iy*lx
            j := (ix+1)%lx+iy*lx

            m.Bond2type[n]   = 0
            m.Bond2hNspin[n] = 2
            m.Bond2weight[n] = 0.5
            m.Tweight += 0.5
            m.Bond2index[n]  = []int{i,j}
            n++
        }
    }
    for iy:=0;iy<ly;iy++ {
        for ix:=0;ix<lx;ix++ {
            i := ix+iy*lx
            j := ix+((iy+1)%ly)*lx

            m.Bond2type[n]   = 0
            m.Bond2hNspin[n] = 2
            m.Bond2weight[n] = 0.5
            m.Tweight += 0.5
            m.Bond2index[n]  = []int{i,j}
            n++
        }
    }
    for iy:=0;iy<ly;iy++ {
        for ix:=0;ix<lx;ix++ {
            i1 := (ix+0)%lx+((iy+0)%ly)*lx
            i2 := (ix+1)%lx+((iy+0)%ly)*lx
            i3 := (ix+0)%lx+((iy+1)%ly)*lx
            i4 := (ix+1)%lx+((iy+1)%ly)*lx
            i5 := (ix+0)%lx+((iy+2)%ly)*lx
            i6 := (ix+1)%lx+((iy+2)%ly)*lx

            m.Bond2type[n]   = 1
            m.Bond2hNspin[n] = 6
            m.Bond2weight[n] = 0.125*q
            m.Tweight += 0.125*q
            m.Bond2index[n]  = []int{i1,i2,i3,i4,i5,i6}
            n++
        }
    }
    for iy:=0;iy<ly;iy++ {
        for ix:=0;ix<lx;ix++ {
            i1 := (ix+0)%lx+((iy+0)%ly)*lx
            i2 := (ix+0)%lx+((iy+1)%ly)*lx
            i3 := (ix+1)%lx+((iy+0)%ly)*lx
            i4 := (ix+1)%lx+((iy+1)%ly)*lx
            i5 := (ix+2)%lx+((iy+0)%ly)*lx
            i6 := (ix+2)%lx+((iy+1)%ly)*lx

            m.Bond2type[n]   = 1
            m.Bond2hNspin[n] = 6
            m.Bond2weight[n] = 0.125*q
            m.Tweight += 0.125*q
            m.Bond2index[n]  = []int{i1,i2,i3,i4,i5,i6}
            n++
        }
    }


    m.LinkRule = make([][]int,2)
    m.LinkRule[0] = link_rule.LinkRuleHeisenbergAFM
    m.LinkRule[1] = link_rule.LinkRuleJQ3

    m.InsertRule = make([]func([]int)bool,2)
    m.InsertRule[0] = insert_rule.InsertRuleHeisenbergAFM
    m.InsertRule[1] = insert_rule.InsertRuleJQ3

    m.Nsite = nsite
    m.Nbond = nbond

    return m
}

func test_models() {
    m := IsingTFSquare(4,4,1.5)

    fmt.Printf("%v\n",m)
}
