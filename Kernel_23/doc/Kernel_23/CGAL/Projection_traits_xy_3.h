
namespace CGAL {

/*!
\ingroup kernel_classes

The class `Projection_traits_xy_3` is an adapter to apply 2D algorithms to the projections of 3D data on the `xy`-plane.

\cgal provides also predefined geometric traits classes
`Projection_traits_yz_3<Gt>` and
`Projection_traits_xz_3<Gt>` to
deal with projections on the
`zx`- and the `zy`-plane,
respectively.

\tparam Gt must be a model of `ProjectionTraitsGeometricTraits_3`

\note Internal constructions (projections) are used in the predicate and
construction functors of this class. If `Gt` is a model of `Kernel` providing exact
constructions or if `Gt` is a `CGAL::Filtered_kernel` (such as for
`CGAL::Exact_predicates_inexact_constructions_kernel`), this class automatically
provides exact predicates.

\cgalModelsBareBegin
\cgalModelsBare{The class is a model of several 2D triangulation traits class concepts,
  except that it does not provide the type and constructors required to build the dual Voronoi diagram.}
\cgalModelsBare{`PolygonTraits_2`}
\cgalModelsBare{`ConvexHullTraits_2`}
\cgalModelsBare{`TriangulationTraits_2`}
\cgalModelsBare{`DelaunayTriangulationTraits_2`}
\cgalModelsBare{`ConstrainedTriangulationTraits_2`}
\cgalModelsBare{`ConvexHullTraits_2`}
\cgalModelsBare{`DelaunayMeshTraits_2`}
\cgalModelsBare{`AnalyticWeightTraits_2`}
\cgalModelsBare{`Barycentric_coordinates::BarycentricTraits_2`}
\cgalModelsBareEnd

\sa `CGAL::Projection_traits_3`
*/
template< typename Gt >
class Projection_traits_xy_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef Point_3<Gt> Point_2;

/*!

*/
typedef Segment_3<Gt> Segment_2;

/*!

*/
typedef Triangle_3<Gt> Triangle_2;

/*!

*/
typedef Line_3<Gt> Line_2;

/// @}

/// \name Functors
/// The functors provided by this class are those listed in the
/// concepts, except that it does not provide the type and
/// constructors required to build the dual Voronoi diagram. The
/// functors operate on the 2D projection of their arguments. They
/// come with preconditions that projections of the arguments are
/// non-degenerate, e.g. a line segment does not project on a single
/// point, two points do not project on the same point, etc. In the
/// following, we specify the choice of the `z`-coordinate in case a
/// new point is constructed.
/// @{

/*!
A construction object.
Provides the operator :

`std::optional< std::variant<Point_2,Segment_2> > operator()(Segment_2 s1, Segment_2 s2);`
which returns a 3D object whose projection on the xy-plane
is the intersection of the projections of `s1` and `s2`.
If non empty, the returned object is either a segment or a point.
Its embedding in 3D is computed as the interpolation
between `s1` and `s2`,
meaning that any point `p` of the returned object
is the midpoint of segment `p1p2` where `p1` and `p2` are the two points of `s1` and `s2` respectively, both projecting on `p`.
\pre The projection of `s1` and the projection of `s2` are non-degenerate `2D` segments.

*/
typedef unspecified_type Intersect_2;

/// @}

/// \name Creation
/// @{

/*!

default constructor.
*/
Projection_traits_xy_3();

/*!
Copy constructor.
*/
Projection_traits_xy_3(
Projection_traits_xy_3 tr);

/*!
Assignment operator.
*/
Projection_traits_xy_3 operator=(Projection_traits_xy_3 tr);

/// @}

}; /* end Projection_traits_xy_3 */
} /* end namespace CGAL */
