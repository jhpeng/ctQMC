package cluster

import (
    "fmt"
    . "github.com/jhpeng/ctQMC/dtype"
)

func Root(p map[Id]Id, v Id) Id {
    u := v
    for p[u]!=u {
        u = p[u]
    }

    return u
}

func Compress(p map[Id]Id, v Id, r Id) {
    u := v
    for p[u]!=u {
        s := p[u]
        p[u] = r
        u = s
    }
}

func Union(p map[Id]Id, w map[Id]uint64, va Id, vb Id){
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
    v0 := Id{Key:0, I:0}
    v1 := Id{Key:1, I:0}
    v2 := Id{Key:2, I:0}
    v3 := Id{Key:3, I:0}
    v4 := Id{Key:4, I:0}
    v5 := Id{Key:5, I:0}
    v6 := Id{Key:6, I:0}
    v7 := Id{Key:7, I:0}

    p := make(map[Id]Id)
    w := make(map[Id]uint64)

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
