package main

import (
    "fmt"
    "time"
    "runtime"
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

var staggeredStructure []float64
func measurement(w WorldLine, m Model, samples []Estimator) {

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
        v       := sequence[i]
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
    lx   = 32
    ly   = 32
    beta = 10.0
    q3   = 1.5
    seed := int64(84686)

    rand.Seed(seed)

    m := models.JQ3LadderSquare(lx,ly,q3)

    w := update.NewWorldLine(12,m.Nsite,beta)

    samples := make([]Estimator,10)

    for i:=0;i<10;i++ {
        samples[i] = stats.New()
    }
    (samples[0]).Name = "ms1"
    (samples[1]).Name = "ms2"
    (samples[2]).Name = "ms4"
    (samples[3]).Name = "Xs"
    (samples[4]).Name = "Xu"
    (samples[5]).Name = "Noo"

    for i:=0;i<w.Nsite;i++ {
        w.State[i] = 1
        if rand.Float64()<0.5 {
            w.State[i] = -1
        }
    }

    t0 := time.Now()
    t1 := time.Now()
    d0 := time.Second*0
    d1 := time.Second*0
    d2 := time.Second*0
    n := 0
    for {
        start := time.Now()
        w = update.Remove(w)
        w = update.Insert(w,m)
        end_graph := time.Now()
        update.InnerLink(w,m)
        update.OuterLink(w,m)
        update.FlipCluster(w)
        end_conf  := time.Now()

        d0 += end_graph.Sub(start)
        d1 += end_conf.Sub(end_graph)

        if n>2000 {
            end_conf = time.Now()
            measurement(w,m,samples)
            t2 := time.Now()
            d2 += t2.Sub(end_conf)

            if t2.Sub(t1)>1*time.Minute {
                t1 = time.Now()
                fmt.Printf("===========================================\n")
                PrintMemUsage()
                fmt.Printf("updating graphs : %v ,  updating configuration : %v ,  measurement : %v\n",d0,d1,d2)
                fmt.Printf("noo=%v nsweep=%v time=%v\n",w.Nvertices,n,t2.Sub(t0))
                stats.PrintDetail(samples[1])
                stats.PrintDetail(samples[3])
                stats.PrintDetail(samples[4])
                stats.PrintDetail(samples[5])
            }
        }

        n++
    }
}
