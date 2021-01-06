package cluster

import (
    "fmt"
)

type Index struct {
    key uint64
    i int
}

func Root(p map[Index]Index, v Index) Index {
    u := v
    for p[u]!=u {
        u = p[u]
    }

    return u
}

func Compress(p map[Index]Index, v Index, r Index) {
    u := v
    for p[u]!=u {
        s := p[u]
        p[u] = r
        u = s
    }
}

func Union(p map[Index]Index, w map[Index]uint64, va Index, vb Index){
    ra := Root(p,va)
    rb := Root(p,vb)
    r0 := ra
    r1 := rb
    if w[ra]<w[rb] {
        r0 = rb
        r1 = ra
    }

    w[r0] = w[r0]+w[r1]
    p[r1] = r0

    Compress(p,va,r0)
    Compress(p,vb,r0)
}

func test() {
    v0 := Index{key:0, i:0}
    v1 := Index{key:1, i:0}
    v2 := Index{key:2, i:0}
    v3 := Index{key:3, i:0}
    v4 := Index{key:4, i:0}
    v5 := Index{key:5, i:0}
    v6 := Index{key:6, i:0}
    v7 := Index{key:7, i:0}

    p := make(map[Index]Index)
    w := make(map[Index]uint64)

    p[v0] = v0
    p[v1] = v1
    p[v2] = v2
    p[v3] = v3
    p[v4] = v4
    p[v5] = v5
    p[v6] = v6
    p[v7] = v7

    w[v0] = 1
    w[v1] = 1
    w[v2] = 1
    w[v3] = 1
    w[v4] = 1
    w[v5] = 1
    w[v6] = 1
    w[v7] = 1

    Union(p,w,v1,v4)
    Union(p,w,v5,v6)
    Union(p,w,v6,v4)

    fmt.Println(p)
    fmt.Println(w)
}
