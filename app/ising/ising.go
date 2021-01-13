package main

import (
    "fmt"
    "time"
    "math/rand"
    "github.com/jhpeng/ctQMC/models"
    "github.com/jhpeng/ctQMC/update"
    . "github.com/jhpeng/ctQMC/dtype"
)

func main() {
    var w WorldLine
    x := 256
    y := 128
    seed := int64(6329012)
    beta := 128.0

    rand.Seed(seed)

    m := models.IsingTFSquare(x,y,1.0)

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
        if rand.Float64()<1.0 {
            w.State[i] = -1
        }
    }

    for i:=0;i<10000;i++ {
        t1 := time.Now()
        w = update.Remove(w)
        w = update.Insert(w,m)
        t2 := time.Now()
        p,c := update.InnerLink(w,m)
        update.OuterLink(w,m,p,c)
        update.FlipCluster(w,p)
        t3 := time.Now()
        fmt.Printf("%v %v %v\n",w.Nvertices,t2.Sub(t1),t3.Sub(t2))
        /*
        for k:=0;k<y;k++ {
            for j:=0;j<x;j++ {
                i_site := j+x*k
                if w.State[i_site]==1 {
                    fmt.Printf("x ")
                } else{
                    fmt.Printf("o ")
                }
            }
            fmt.Printf("\n")
        }*/
    }
}
