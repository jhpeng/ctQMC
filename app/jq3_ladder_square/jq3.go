package main

import (
    "fmt"
    "time"
    //"runtime"
    "math/rand"
    "github.com/jhpeng/ctQMC/models"
    "github.com/jhpeng/ctQMC/update"
    . "github.com/jhpeng/ctQMC/dtype"
)

func numberOfGraphs(w WorldLine, m Model) {
    var key uint64

    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    n0:=0
    n1:=0

    for i:=0;i<w.Nvertices;i++ {
        key = sequence[i]
        v    := w.Table[key]
        bond := v.Bond
        t    := m.Bond2type[bond]
        if t==0 {
            n0++
        } else if t==1 {
            n1++
        } else {
            fmt.Printf("something wrong!\n")
        }
    }

    fmt.Printf("n0=%v n1=%v\n",n0,n1)
}

func checkPeriodic(w WorldLine, m Model) bool {
    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    pstate := make([]int,w.Nsite)
    for i:=0;i<w.Nsite;i++ {
        pstate[i] = w.State[i]
    }

    for i:=0;i<w.Nvertices;i++ {
        key     := sequence[i]
        v       := w.Table[key]
        bond    := v.Bond
        hNspin  := v.HNspin
        indices := m.Bond2index[bond]

        for j:=0;j<hNspin;j++ {
            i_site := indices[j]
            if pstate[i_site]!=v.State[j]{
                fmt.Printf("something wrong!\n")
            }
            pstate[i_site] = v.State[hNspin+j]
        }
    }

    check := true
    for i:=0;i<w.Nsite;i++ {
        if pstate[i]!=w.State[i] {
            check = false
        }
    }

    return check
}

func main() {
    var w WorldLine
    x := 32
    y := 32
    beta := 64.0
    q := 1.5
    seed := int64(32901)

    rand.Seed(seed)

    m := models.JQ3LadderSquare(x,y,q)
    //fmt.Println(m)

    w.Table = make(map[uint64]Vertex)
    w.SequenceA = make([]uint64,2048)
    w.SequenceB = make([]uint64,2048)

    w.Nvertices = 0
    w.Flag  = true
    w.State = make([]int,m.Nsite)
    w.Last  = make([]Id,m.Nsite)
    w.First = make([]Id,m.Nsite)
    w.Nsite = m.Nsite
    w.Beta = beta

    for i:=0;i<w.Nsite;i++ {
        w.State[i] = 1
        if rand.Float64()<0.5 {
            w.State[i] = -1
        }
    }

    t := time.Now()
    n := 0
    for {
        t1  := time.Now()
        w    = update.Remove(w)
        w    = update.Insert(w,m)
        t2  := time.Now()
        p,c := update.InnerLink(w,m)
        update.OuterLink(w,m,p,c)
        update.FlipCluster(w,p)
        t3 := time.Now()
        if t3.Sub(t)>time.Second {
            t = time.Now()
            fmt.Printf("noo=%v %v %v nsweep=%v\n",w.Nvertices,t2.Sub(t1),t3.Sub(t2),n)
            //numberOfGraphs(w,m)
        }
        //if !checkPeriodic(w,m) {
        //    fmt.Printf("false!\n")
        //}
        n++
    }
}
