use parry3d_f64::{
    math::{Pose, Rotation, Vector},
    query::{contact, intersection_test, PointQuery},
    shape::{Cuboid, Cylinder, SupportMap, Triangle},
};

#[test]
fn coincident_nested_cylinders_intersect_in_both_orders() {
    let outer = Cylinder::new(50.0, 100.0);
    let inner = Cylinder::new(1.5, 1.0);
    assert!(intersection_test(&Pose::IDENTITY, &outer, &Pose::IDENTITY, &inner).unwrap());
    assert!(intersection_test(&Pose::IDENTITY, &inner, &Pose::IDENTITY, &outer).unwrap());
    let separated = Pose::from_translation(Vector::new(102.0, 0.0, 0.0));
    assert!(!intersection_test(&Pose::IDENTITY, &outer, &separated, &inner).unwrap());
}

#[test]
fn large_triangle_cylinder_edge_contact_retains_native_witnesses() {
    use parry3d_f64::query::{
        details::contact_support_map_support_map_with_params,
        gjk::{GJKResult, VoronoiSimplex},
    };
    let triangle = Triangle::new(
        Vector::new(-250.0, 0.0, -250.0),
        Vector::new(250.0, 0.0, 250.0),
        Vector::new(250.0, 0.0, -250.0),
    );
    let cylinder = Cylinder::new(0.15, 0.4);
    let pose = Pose::from_parts(
        Vector::new(26.275517120869555, 0.3989396412921885, 26.459245381148524),
        Rotation::from_xyzw(
            0.10039221371651696,
            0.1742794765226331,
            0.6997061011501934,
            -0.6855358526530525,
        ),
    );
    let axis = Vector::new(
        -0.09186413013952333,
        0.3989396412921885,
        0.09186413013944517,
    )
    .normalize();
    let result = contact_support_map_support_map_with_params(
        &pose,
        &triangle,
        &cylinder,
        496.8394434529555,
        &mut VoronoiSimplex::new(),
        Some(axis),
    );
    let GJKResult::ClosestPoints(first, second, normal) = result else {
        panic!("enclosing prediction lost native geometry: {result:?}")
    };
    assert!(first.is_finite() && second.is_finite() && normal.is_finite());
    assert!((normal.length() - 1.0).abs() < 1e-12);
    assert!((triangle.project_local_point(first, false).point - first).length() < 1e-9);
    let boundary_error = (cylinder.project_point(&pose, second, false).point - second).length();
    // GJK witnesses interpolate supports on a curved boundary. Their surface
    // residual follows its relative convergence precision, not exact vertex arithmetic.
    let tolerance = 4.0 * parry3d_f64::query::gjk::eps_tol().sqrt() * cylinder.radius;
    assert!(
        boundary_error < tolerance,
        "cylinder boundary error {boundary_error:e}"
    );
    let distance = (second - first).dot(normal);
    assert!(distance > 0.0 && distance < 0.002);
    assert!((second - first - normal * distance).length() < tolerance);
}

#[test]
fn cylinder_interior_projection_selects_a_nearest_boundary_on_cap_ties() {
    let cylinder = Cylinder::new(0.1, 0.2);
    for point in [Vector::ZERO, Vector::new(0.05, 0.0, 0.0)] {
        let projected = cylinder.project_local_point(point, false);
        assert!(projected.is_inside);
        assert!(((projected.point - point).length() - 0.1).abs() < 1e-12);
        assert!((projected.point.y.abs() - 0.1).abs() < 1e-12);
        assert_eq!(cylinder.project_local_point(point, true).point, point);
    }
    let tall = Cylinder::new(1.0, 0.2);
    assert!((tall.project_local_point(Vector::ZERO, false).point.length() - 0.2).abs() < 1e-12);
}

#[test]
fn touching_triangle_cuboid_witnesses_belong_to_their_shapes() {
    let vertex = Vector::new(0.0, -2.0, -2.0);
    let triangle = Triangle::new(vertex + Vector::Z * 4.0, vertex + Vector::Y * 4.0, vertex);
    let cuboid = Cuboid::new(Vector::new(0.125, 0.25, 0.25));
    let normal = Vector::new(1.0, -1.0, -1.0).normalize();
    let pose = Pose::from_translation(vertex - cuboid.local_support_point(-normal));
    let hit = contact(&Pose::IDENTITY, &triangle, &pose, &cuboid, 1.0)
        .unwrap()
        .unwrap();
    assert!(hit.dist.abs() < 1e-12);
    assert!((hit.point1 - vertex).length() < 1e-12);
    assert!((hit.point2 - vertex).length() < 1e-12);
    assert!((hit.point2 - hit.point1 - hit.normal1 * hit.dist).length() < 1e-12);
    let reversed = contact(&pose, &cuboid, &Pose::IDENTITY, &triangle, 1.0)
        .unwrap()
        .unwrap();
    assert!((reversed.point1 - hit.point2).length() < 1e-12);
    assert!((reversed.point2 - hit.point1).length() < 1e-12);
    let gap = Pose::from_translation(pose.translation + normal * 0.01);
    assert!(contact(&Pose::IDENTITY, &triangle, &gap, &cuboid, 0.001)
        .unwrap()
        .is_none());
}
