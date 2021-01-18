package main

import (
    "fmt"
    "time"
    //"runtime"
    "math"
    "math/rand"
    "github.com/jhpeng/ctQMC/models"
    "github.com/jhpeng/ctQMC/update"
    "github.com/jhpeng/ctQMC/stats"
    . "github.com/jhpeng/ctQMC/dtype"
)

var (
    lx int
    ly int
    beta float64
    q3 float64
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

var staggeredStructure []float64
func measurement(w WorldLine, m Model, samples []Estimator) {
    var key uint64

    if staggeredStructure==nil {
        staggeredStructure = make([]float64,w.Nsite)
        for j:=0;j<ly;j++ {
            for i:=0;i<lx;i++ {
                i_site := i+j*lx
                staggeredStructure[i_site] = float64(((i+j)%2)*2-1)
            }
        }
    }

    var ( 
        mz float64 = 0
        ms float64 = 0
        msx float64 = 0
        ms1 float64 = 0
        ms2 float64 = 0
        ms4 float64 = 0
        chiu float64 = 0
    )

    for i:=0;i<w.Nsite;i++ {
        mz += float64(w.State[i])
        ms += float64(w.State[i])*staggeredStructure[i]
    }

    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    for i:=0;i<w.Nvertices;i++ {
        key      = sequence[i]
        v       := w.Table[key]
        bond    := v.Bond
        hNspin  := v.HNspin
        state   := v.State
        indices := m.Bond2index[bond]

        for j:=0;j<hNspin;j++ {
            i_site := indices[j]
            dif    := float64((state[j]*state[j+hNspin]-1)*state[j])
            ms += staggeredStructure[i_site]*dif
        }
        msx += ms
        ms1 += math.Abs(ms)
        ms2 += ms*ms
        ms4 += ms*ms*ms*ms
    }

    l   := float64(w.Nvertices)
    vol := float64(w.Nsite)

    chiu = mz*mz/float64(w.Nsite)*beta*0.25
    if w.Nvertices!=0 {
        msx = beta*(msx*msx+ms2)/l/(l+1)/vol*0.25
        ms1 = ms1*0.5/vol/l
        ms2 = ms2*0.25/vol/vol/l
        ms4 = ms4*0.0625/vol/vol/vol/vol/l
    } else {
        msx = beta*ms*ms/2.0/vol*0.25
        ms1 = math.Abs(ms)*0.5/vol
        ms2 = ms*ms*0.25/vol/vol
        ms4 = ms*ms*ms*ms/0.0625/vol/vol/vol/vol
    }

    samples[0] = stats.Append(samples[0],ms1)
    samples[1] = stats.Append(samples[1],ms2)
    samples[2] = stats.Append(samples[2],ms4)
    samples[3] = stats.Append(samples[3],msx)
    samples[4] = stats.Append(samples[4],chiu)
    samples[5] = stats.Append(samples[5],float64(w.Nvertices))
}

func main() {
    var w WorldLine
    lx   = 32
    ly   = 32
    beta = 10.0
    q3   = 1.5
    seed := int64(783621)

    rand.Seed(seed)

    m := models.JQ3LadderSquare(lx,ly,q3)
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

    samples := make([]Estimator,10)

    for i:=0;i<10;i++ {
        samples[i] = stats.New()
    }

    for i:=0;i<w.Nsite;i++ {
        w.State[i] = 1
        if rand.Float64()<0.5 {
            w.State[i] = -1
        }
    }

    t := time.Now()
    n := 0
    for {
        w    = update.Remove(w)
        w    = update.Insert(w,m)
        p,c := update.InnerLink(w,m)
        update.OuterLink(w,m,p,c)
        update.FlipCluster(w,p)

        if n>1000 {
            measurement(w,m,samples)
        }

        if (n>1000) && (n%1000==0) {
            t1 := time.Now()
            fmt.Printf("--------------------------------------\n")
            fmt.Printf("noo=%v nsweep=%v time=%v\n",w.Nvertices,n,t1.Sub(t))
            fmt.Printf("ms2  = %v\n",stats.Mean(samples[1]))
            fmt.Printf("chis = %v\n",stats.Mean(samples[3]))
            fmt.Printf("chiu = %v\n",stats.Mean(samples[4]))
            fmt.Printf("noo  = %v\n",stats.Mean(samples[5]))
            t = t1
        }

        n++
    }
}
