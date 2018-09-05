// Copyright (c) 2011 Mateusz Czapliński (Go port)
// Copyright (c) 2011 Mahir Iqbal (as3 version)
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// based on http://code.google.com/p/as3polyclip/ (MIT licensed)
// and code by Martínez et al: http://wwwdi.ujaen.es/~fmartin/bool_op.html (public domain)

// Package polyclip provides implementation of algorithm for Boolean operations on 2D polygons.
// For further details, consult the description of Polygon.Construct method.
package polyclip

import (
	"encoding/json"
	"fmt"
	"math"
)

type Point struct {
	X, Y float64
}

func (p *Point) MarshalJSON() ([]byte, error) {
	return []byte(`[` + fmt.Sprintf("%.8f", p.X) + `,` + fmt.Sprintf("%.8f", p.Y) + `]`), nil
}

type point Point

func (p *Point) UnmarshalJSON(b []byte) (err error) {
	var j [2]float64 // := make(map[string]float64)
	err = json.Unmarshal(b, &j)
	if err != nil { //Point é um map
		j := make(map[string]float64)
		err = json.Unmarshal(b, &j)
		if err != nil {
			panic(err)
		}
		p.X = j["x"]
		p.Y = j["y"]
		return
	}
	if err == nil {
		p.X = j[0]
		p.Y = j[1]
		return
	}
	panic(err)
	return
}

// Equals returns true if both p1 and p2 describe exactly the same point.
func (p1 Point) Equals(p2 Point) bool {
	return p1.X == p2.X && p1.Y == p2.Y
}

// Length returns distance from p to point (0, 0).
func (p Point) Length() float64 {
	return math.Sqrt(p.X*p.X + p.Y*p.Y)
}

type Rectangle struct {
	Min, Max Point
}

func (r1 Rectangle) union(r2 Rectangle) Rectangle {
	return Rectangle{
		Min: Point{
			X: math.Min(r1.Min.X, r2.Min.X),
			Y: math.Min(r1.Min.Y, r2.Min.Y),
		},
		Max: Point{
			X: math.Max(r1.Max.X, r2.Max.X),
			Y: math.Max(r1.Max.Y, r2.Max.Y),
		}}
}

// Overlaps returns whether r1 and r2 have a non-empty intersection.
func (r1 Rectangle) Overlaps(r2 Rectangle) bool {
	return r1.Min.X <= r2.Max.X && r1.Max.X >= r2.Min.X &&
		r1.Min.Y <= r2.Max.Y && r1.Max.Y >= r2.Min.Y
}

// Return the rectangle's center
func (r1 Rectangle) Center() Point {
	return Point{(r1.Max.X + r1.Min.X) / 2.0, (r1.Max.Y + r1.Min.Y) / 2.0}
}

// Used to represent an edge of a polygon.
type segment struct {
	start, end Point
}

// Contour represents a sequence of vertices connected by line segments, forming a closed shape.
type Contour []Point

// Add is a convenience method for appending a point to a contour.
func (c *Contour) Add(p Point) {
	*c = append(*c, p)
}

// BoundingBox finds minimum and maximum coordinates of points in a contour.
func (c Contour) BoundingBox() Rectangle {
	bb := Rectangle{}
	bb.Min.X = math.Inf(1)
	bb.Min.Y = math.Inf(1)
	bb.Max.X = math.Inf(-1)
	bb.Max.Y = math.Inf(-1)

	for _, p := range c {
		if p.X > bb.Max.X {
			bb.Max.X = p.X
		}
		if p.X < bb.Min.X {
			bb.Min.X = p.X
		}
		if p.Y > bb.Max.Y {
			bb.Max.Y = p.Y
		}
		if p.Y < bb.Min.Y {
			bb.Min.Y = p.Y
		}
	}
	return bb
}

func (c Contour) segment(index int) segment {
	if index == len(c)-1 {
		return segment{c[len(c)-1], c[0]}
	}
	return segment{c[index], c[index+1]}
	// if out-of-bounds, we expect panic detected by runtime
}

// Checks if a point is inside a contour using the "point in polygon" raycast method.
// This works for all polygons, whether they are clockwise or counter clockwise,
// convex or concave.
// See: http://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
// Returns true if p is inside the polygon defined by contour.
func (c Contour) Contains(p Point) bool {
	// Cast ray from p.x towards the right
	intersections := 0
	for i := range c {
		curr := c[i]
		ii := i + 1
		if ii == len(c) {
			ii = 0
		}
		next := c[ii]

		if (p.Y >= next.Y || p.Y <= curr.Y) &&
			(p.Y >= curr.Y || p.Y <= next.Y) {
			continue
		}
		// Edge is from curr to next.

		if p.X >= math.Max(curr.X, next.X) ||
			next.Y == curr.Y {
			continue
		}

		// Find where the line intersects...
		xint := (p.Y-curr.Y)*(next.X-curr.X)/(next.Y-curr.Y) + curr.X
		if curr.X != next.X && p.X > xint {
			continue
		}

		intersections++
	}

	return intersections%2 != 0
}

// Checks if a point is inside a contour using the "point in polygon" raycast method.
// This works for all polygons, whether they are clockwise or counter clockwise,
// convex or concave.
// See: http://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
// Returns true if p is inside the polygon defined by contour.
func (c Contour) ContainsWnPoly(p Point) bool {
	cn := 0
	for i := range c { // edge from c[i]  to nextC
		C := c[i]
		var nextC Point
		if i+1 == len(c) {
			nextC = c[0]
		} else {
			nextC = c[i+1]
		}
		if ((C.Y <= p.Y) && (nextC.Y > p.Y)) || ((C.Y > p.Y) && (nextC.Y <= p.Y)) { // a downward crossing
			// compute  the actual edge-ray intersect x-coordinate
			vt := float64(p.Y-C.Y) / (nextC.Y - C.Y)
			if p.X < C.X+vt*(nextC.X-C.X) { // p.Coordinates.X < intersect
				cn++ // a valid crossing of y=p.Coordinates[1] right of p.Coordinates[0]
			}
		}
	}
	return (cn&1 == 1) // 0 if even (out), and 1 if  odd (in)
}

func (c Contour) Smooth() {
	prevPt := c[0]
	nPoints := len(c)
	for i, pt := range c[1:] { // edge from c[i]  to nextC
		pt.X += prevPt.X
		pt.X /= 2.0
		pt.Y += prevPt.Y
		pt.Y /= 2.0
		prevPt = pt
		c[i] = pt
	}
	c[nPoints-1] = c[0]
}

func (c Contour) GetReversed() Contour {
	var rc Contour
	L := len(c) - 1
	for i, _ := range c { // edge from c[i]  to nextC
		rc.Add(c[L-i])
	}
	return rc
}

func (c Contour) SmoothRDP(epsilon float64) Contour {
	a := smoothRDP(c, epsilon)
	a[len(a)-1] = a[0]
	return a
}

func rDPfindPerpendicularDistance(p, p1, p2 Point) (result float64) {
	if p1.X == p2.X {
		result = math.Abs(float64(p.X - p1.X))
	} else {
		slope := float64(p2.Y-p1.Y) / float64(p2.X-p1.X)
		intercept := float64(p1.Y) - (slope * float64(p1.X))
		result = math.Abs(slope*float64(p.X)-float64(p.Y)+intercept) / math.Sqrt(math.Pow(slope, 2)+1)
	}
	return
}

func smoothRDP(points []Point, epsilon float64) []Point {
	if len(points) < 3 {
		return points
	}
	firstPoint := points[0]
	lastPoint := points[len(points)-2]
	index := -1
	dist := float64(0)
	for i := 1; i < len(points)-1; i++ {
		cDist := rDPfindPerpendicularDistance(points[i], firstPoint, lastPoint)
		if cDist > dist {
			dist = cDist
			index = i
		}
	}
	if dist > epsilon {
		l1 := points[0 : index+1]
		l2 := points[index:]
		r1 := smoothRDP(l1, epsilon)
		r2 := smoothRDP(l2, epsilon)
		rs := append(r1[0:len(r1)-1], r2...)
		return rs
	} else {
		ret := make([]Point, 0)
		ret = append(ret, firstPoint, lastPoint)
		return ret
	}
	return make([]Point, 0)
}

func (c Contour) SignedArea() float64 {
	sA := 0.0
	prev := c[0]
	for _, pt := range c[1:] { // edge from c[i]  to nextC
		sA += (prev.X*pt.Y - prev.Y*pt.X)
		prev = pt
	}
	return sA / 2.0
}

// Clone returns a copy of a contour.
func (c Contour) Clone() Contour {
	return append([]Point{}, c...)
}

// Polygon is carved out of a 2D plane by a set of (possibly disjoint) contours.
// It can thus contain holes, and can be self-intersecting.
type Polygon []Contour

// NumVertices returns total number of all vertices of all contours of a polygon.
func (p Polygon) NumVertices() int {
	num := 0
	for _, c := range p {
		num += len(c)
	}
	return num
}

// BoundingBox finds minimum and maximum coordinates of points in a polygon.
func (p Polygon) BoundingBox() Rectangle {
	//fmt.Printf("\nPolygon pp->%#v\n", p)
	bb := p[0].BoundingBox()
	for _, c := range p[1:] {
		bb = bb.union(c.BoundingBox())
	}

	return bb
}

// Add is a convenience method for appending a contour to a polygon.
func (p *Polygon) Add(c Contour) {
	*p = append(*p, c)
}

// Clone returns a duplicate of a polygon.
func (p Polygon) Clone() Polygon {
	r := Polygon(make([]Contour, len(p)))
	for i := range p {
		r[i] = p[i].Clone()
	}
	return r
}

// Op describes an operation which can be performed on two polygons.
type Op int

const (
	UNION Op = iota
	INTERSECTION
	DIFFERENCE
	XOR
)

// Construct computes a 2D polygon, which is a result of performing
// specified Boolean operation on the provided pair of polygons (p <Op> clipping).
// It uses algorithm described by F. Martínez, A. J. Rueda, F. R. Feito
// in "A new algorithm for computing Boolean operations on polygons"
// - see: http://wwwdi.ujaen.es/~fmartin/bool_op.html
// The paper describes the algorithm as performing in time O((n+k) log n),
// where n is number of all edges of all polygons in operation, and
// k is number of intersections of all polygon edges.
func (p Polygon) Construct(operation Op, clipping Polygon) Polygon {
	c := clipper{
		subject:  p,
		clipping: clipping,
	}
	return c.compute(operation)
}
