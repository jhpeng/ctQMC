package main

import (
    "fmt"
    "time"
    "runtime"
    "math/rand"
    "github.com/jhpeng/ctQMC/models"
    "github.com/jhpeng/ctQMC/update"
    . "github.com/jhpeng/ctQMC/dtype"
)

func bToMb(b uint64) uint64 {
    return b / 1024 / 1024
}

func PrintMemUsage() {
    var m runtime.MemStats
    runtime.ReadMemStats(&m)
    // For info on each, see: https://golang.org/pkg/runtime/#MemStats
    fmt.Printf("Alloc = %v MiB", bToMb(m.Alloc))
    fmt.Printf("\tTotalAlloc = %v MiB", bToMb(m.TotalAlloc))
    fmt.Printf("\tSys = %v MiB", bToMb(m.Sys))
    fmt.Printf("\tNumGC = %v\n", m.NumGC)
}

func measurement(w WorldLine) {
    var mz float64 = 0.0

    for i:=0;i<w.Nsite;i++ {
        mz += float64(w.State[i])
    }

    mz = mz/float64(w.Nsite)
    mz2 := mz*mz

    fmt.Printf("mz  = %v \n",mz)
    fmt.Printf("mz2 = %v \n",mz2)
}

func main() {
    var w WorldLine
    x := 1024
    y := 1024
    seed := int64(329012)
    beta := 1.80

    rand.Seed(seed)

    m := models.IsingTFSquare(x,y,0.5)

    length := 2048
    w.SequenceA = make([]Vertex,length)
    w.SequenceB = make([]Vertex,length)
    w.Mnspin    = 4
    w.Cluster   = make([]int,w.Mnspin*length)
    w.Weight    = make([]int,w.Mnspin*length)

    w.Nvertices = 0
    w.Flag  = true
    w.State = make([]int,m.Nsite)
    w.Last  = make([]int,m.Nsite)
    w.First = make([]int,m.Nsite)
    w.Nsite = m.Nsite
    w.Beta = beta

    for i:=0;i<w.Nsite;i++ {
        w.State[i] = 1
        if rand.Float64()<0.5 {
            w.State[i] = -1
        }
    }

    for {
        t1 := time.Now()
        w = update.Remove(w)
        w = update.Insert(w,m)
        t2 := time.Now()
        update.InnerLink(w,m)
        update.OuterLink(w,m)
        update.FlipCluster(w)
        t3 := time.Now()
        fmt.Printf("%v %v %v\n",w.Nvertices,t2.Sub(t1),t3.Sub(t2))
        measurement(w)
        PrintMemUsage()

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
