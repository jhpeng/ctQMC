package update

import (
    "fmt"
    "math"
    "math/rand"
    "github.com/jhpeng/ctQMC/update/cluster"
    . "github.com/jhpeng/ctQMC/dtype"
)

var cmf []float64
var count []uint64
func pointSamplingWeighted(t float64, w []float64) int {
    dis := rand.Float64()*t

    if cmf==nil {
        cmf    = make([]float64, len(w))
        count  = make([]uint64,  len(w))
        cmf[0] = w[0]
        count[0] = 0
        for i:=1;i<len(w);i++ {
            cmf[i] = cmf[i-1]+w[i]
            count[i] = 0
        }
    }

    i:=-1
    j:=len(cmf)
    d:=j-i
    for d>1 {
        if dis > cmf[i+d/2] {
            i = i+d/2
        } else {
            j = i+d/2
        }
        d = j-i
    }

    count[j]++

    return j
}

var MaxKeyFromSample uint64 = 1000000000000000 // 10E15
func squenceSamplingUniform(lambda float64, start uint64) []uint64{
    var sequence []uint64
    var k float64 = 0

    dis := rand.Float64()
    k -= math.Log(dis)/lambda
    for k<1.0 {
        sequence = append(sequence,uint64(k*float64(MaxKeyFromSample))+start)
        dis = rand.Float64()
        k -= math.Log(dis)/lambda
    }

    return sequence
}

func NewWorldLine(mnspin int, nsite int, beta float64) WorldLine {
    var w WorldLine

    length := 1024
    w.SequenceA = make([]Vertex,length)
    w.SequenceB = make([]Vertex,length)
    w.Mnspin    = mnspin
    w.Cluster   = make([]int,mnspin*length)
    w.Weight    = make([]int,mnspin*length)

    w.Nvertices = 0
    w.Flag  = true
    w.State = make([]int,nsite)
    w.Last  = make([]int,nsite)
    w.First = make([]int,nsite)
    w.Nsite = nsite
    w.Beta = beta

    return w
}

func Remove(w WorldLine) WorldLine{
    var v   Vertex

    sequence1 := w.SequenceB
    sequence2 := w.SequenceA
    if w.Flag {
        sequence1 = w.SequenceA
        sequence2 = w.SequenceB
    }

    k:=0
    for i:=0;i<w.Nvertices;i++ {
        v = sequence1[i]
        check_delete := true

        for j:=0;j<v.HNspin;j++ {
            if v.State[j]!=v.State[j+v.HNspin]{
                check_delete = false
            }
        }

        if !check_delete {
            sequence2[k] = sequence1[i]
            k++
        }
    }
    w.Nvertices = k
    w.Flag = !(w.Flag)

    return w
}

func Insert(w WorldLine, m Model) WorldLine{
    lambda := (m.Tweight)*(w.Beta)
    insert_seq := squenceSamplingUniform(lambda,0)

    length := len(insert_seq)+w.Nvertices+1024
    if length>len(w.SequenceA) {
        temp1 := make([]Vertex,length)
        copy(temp1,w.SequenceA)
        temp2 := make([]Vertex,length)
        copy(temp2,w.SequenceB)
        w.SequenceA = temp1
        w.SequenceB = temp2

        w.Cluster = make([]int,length*w.Mnspin)
        w.Weight  = make([]int,length*w.Mnspin)
    }

    sequence1 := w.SequenceB
    sequence2 := w.SequenceA
    if w.Flag {
        sequence1 = w.SequenceA
        sequence2 = w.SequenceB
    }

    pstate := make([]int,m.Nsite)
    copy(pstate,w.State)

    k:=0
    n:=0
    var v Vertex
    for i:=0;i<len(insert_seq);i++ {
        key1 := (sequence1[k]).Key
        key2 := insert_seq[i]

        for (key1<key2) && (k<w.Nvertices) {
            v = sequence1[k]
            indices := m.Bond2index[v.Bond]
            for i_site:=0;i_site<v.HNspin;i_site++ {
                index := indices[i_site]
                pstate[index] = v.State[v.HNspin+i_site]
            }

            sequence2[n] = v
            n++
            k++
            key1 = (sequence1[k]).Key
        }

        if key1!=key2 {
            bond    := pointSamplingWeighted(m.Tweight,m.Bond2weight)
            t       := m.Bond2type[bond]
            indices := m.Bond2index[bond]
            hNspin  := m.Bond2hNspin[bond]
            rule    := m.InsertRule[t]
            lstate  := make([]int,2*hNspin)
            for i_site:=0;i_site<hNspin;i_site++ {
                lstate[i_site]        = pstate[indices[i_site]]
                lstate[i_site+hNspin] = pstate[indices[i_site]]
            }

            if rule(lstate) {
                var iv Vertex
                iv.Key    = key2
                iv.Bond   = bond
                iv.HNspin = hNspin
                iv.State  = lstate

                sequence2[n] = iv
                n++
            }
        }
    }

    for k<w.Nvertices {
        sequence2[n] = sequence1[k]
        n++
        k++
    }

    w.Nvertices = n
    w.Flag = !(w.Flag)

    return w
}

func InnerLink(w WorldLine, m Model) {

    var (
        v      Vertex
        bond   int
        hNspin int
        t      int
        idn    int
        idv    int
    )

    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    for i:=0;i<w.Nvertices;i++ {
        v      = sequence[i]
        bond   = v.Bond
        hNspin = v.HNspin
        t      = m.Bond2type[bond]
        rule  := m.LinkRule[t]

        for j:=0;j<2*hNspin;j++ {
            idn = i*w.Mnspin+j
            idv = i*w.Mnspin+rule[j]
            w.Cluster[idn] = idv
            w.Weight[idn]  = rule[j+2*hNspin]
        }
    }
}

func OuterLink(w WorldLine, m Model) {

    var (
        v         Vertex
        bond      int
        hNspin    int
        indices []int
        index     int
        idp       int
        idn       int
    )

    last  := w.Last
    first := w.First

    for i:=0;i<w.Nsite;i++ {
        last[i]  = -1
        first[i] = -1
    }

    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    for i:=0;i<w.Nvertices;i++ {
        v       = sequence[i]
        bond    = v.Bond
        hNspin  = v.HNspin
        indices = m.Bond2index[bond]

        for j:=0;j<hNspin;j++ {
            index = indices[j]
            idp = i*w.Mnspin+j
            idn = i*w.Mnspin+j+hNspin
            if first[index]==-1 {
                first[index] = idp
                last[index]  = idn
            } else {
                cluster.Union(w.Cluster,w.Weight,last[index],idp)
                last[index] = idn
            }
        }
    }

    for i:=0;i<w.Nsite;i++ {
        if first[i]!=-1 {
            cluster.Union(w.Cluster,w.Weight,first[i],last[i])
        }
    }
}

func FlipCluster(w WorldLine) {

    var(
        v       Vertex
        hNspin  int
        state []int
        idv     int
        idr     int
        id      int
        p       int
        i       int
        j       int
    )

    newState := make(map[int]int)
    first    := w.First

    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    for i=0;i<w.Nvertices;i++ {
        v       = sequence[i]
        hNspin  = v.HNspin
        state   = v.State

        for j=0;j<2*hNspin;j++ {
            idv = i*w.Mnspin+j
            idr = cluster.Root(w.Cluster,idv)
            if newState[idr]==0 {
                if rand.Float64()<0.5 {
                    newState[idr]= 1
                } else {
                    newState[idr]=-1
                }
            }

            state[j] = state[j]*newState[idr]
        }
    }

    for i=0;i<w.Nsite;i++ {
        id  = first[i]
        if id!=-1 {
            p = id/w.Mnspin
            j = id%w.Mnspin

            v = sequence[p]
            w.State[i] = v.State[j]
        } else if rand.Float64()<0.5 {
            w.State[i] =  1
        } else {
            w.State[i] = -1
        }
    }
}

func CheckPeriodic(w WorldLine, m Model) bool {
    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    pstate := make([]int,w.Nsite)
    for i:=0;i<w.Nsite;i++ {
        pstate[i] = w.State[i]
    }

    for i:=0;i<w.Nvertices;i++ {
        v       := sequence[i]
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

