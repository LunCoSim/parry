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

fn check_triangle_cylinder(name: &str, triangle: Triangle, cylinder: Cylinder, pose: Pose) {
    use parry3d_f64::query::{DefaultQueryDispatcher, QueryDispatcher};
    use parry3d_f64::shape::Shape;
    // Equivalent shape-local translations must not change contact existence.
    for origin in [
        Vector::ZERO,
        triangle.center(),
        triangle.project_local_point(pose.translation, false).point,
    ] {
        let triangle = triangle.transformed(&Pose::from_translation(-origin));
        let pose = Pose::from_parts(pose.translation - origin, pose.rotation.normalize());
        for swapped in [false, true] {
            let (a, b, relative): (&dyn Shape, &dyn Shape, Pose) = if swapped {
                (&cylinder, &triangle, pose.inverse())
            } else {
                (&triangle, &cylinder, pose)
            };
            let contact = DefaultQueryDispatcher
                .contact(&relative, a, b, f64::MAX)
                .unwrap()
                .unwrap_or_else(|| panic!("lost {name} origin={origin:?} swapped={swapped}"));
            assert!(
                contact.dist > 0.0 && contact.dist < 0.0001,
                "{name}: {contact:?}"
            );
            assert!((contact.normal1.length() - 1.0).abs() < 1e-12);
            assert!(
                (a.project_local_point(contact.point1, false).point - contact.point1).length()
                    < 1e-6
            );
            assert!(
                (b.project_local_point(contact.point2, false).point - contact.point2).length()
                    < 1e-6
            );
            let separation = relative * contact.point2 - contact.point1;
            assert!((separation.dot(contact.normal1) - contact.dist).abs() < 1e-12);
            assert!((separation - contact.normal1 * contact.dist).length() < 1e-6);
            assert!(DefaultQueryDispatcher
                .contact(&relative, a, b, contact.dist * 0.5)
                .unwrap()
                .is_none());
        }
    }
}

#[test]
fn thin_triangle_cylinder_simplexes_do_not_invent_intersections() {
    use parry3d_f64::math::{Rotation as R, Vector as V};
    check_triangle_cylinder(
        "new-ramp",
        Triangle::new(
            V::new(-7., 1.5941457748413086, -84.245361328125),
            V::new(7., 1.6570005416870117, -84.43522644042969),
            V::new(-7., 1.6570005416870117, -84.43522644042969),
        ),
        Cylinder::new(0.15, 0.34),
        Pose::from_parts(
            V::new(-2.901950571207935, 1.933080902938535, -84.18723223276312),
            R::from_xyzw(
                -0.6959158512406957,
                -0.7030484940800048,
                -0.08118965418870049,
                0.12178730963424969,
            ),
        ),
    );
    check_triangle_cylinder(
        "ramp",
        Triangle::new(
            V::new(-7., 27.257265090942383, -100.2407455444336),
            V::new(7., 27.45724105834961, -100.2437973022461),
            V::new(-7., 27.45724105834961, -100.2437973022461),
        ),
        Cylinder::new(0.125, 0.34),
        Pose::from_parts(
            V::new(-0.4970461203727527, 27.323663729637044, -99.90212374436804),
            R::from_xyzw(
                0.24569900216185098,
                0.23935767672281133,
                0.6628784353589657,
                -0.6655314031944974,
            ),
        ),
    );
    check_triangle_cylinder(
        "offset",
        Triangle::new(
            V::new(14., -1.22340751442872, -14.),
            V::new(14.5, -1.2294075144287198, -14.),
            V::new(14.5, -1.26440751442872, -14.5),
        ),
        Cylinder::new(0.02, 0.34999997913837433),
        Pose::from_parts(
            V::new(14.104314712558768, -1.2037191057296206, -14.083198354068372),
            R::from_xyzw(
                -0.026920021362250367,
                -0.9226774210591173,
                0.03544966899621409,
                0.38299479119906205,
            ),
        ),
    );
}

#[test]
fn thin_tetrahedron_separating_face_survives_simplex_reduction() {
    use parry3d_f64::query::gjk::{CsoPoint, VoronoiSimplex};
    let vertices = [
        Vector::new(
            0.16086747602049556,
            -0.013191043519555468,
            -0.1608630182541333,
        ),
        Vector::new(
            -0.33914823770131086,
            0.02780808043638139,
            0.33912272646340824,
        ),
        Vector::new(
            0.16087233407619017,
            -0.013190772712001975,
            -0.16085861144057745,
        ),
        Vector::new(
            0.1608708775151637,
            -0.013190853905066247,
            -0.16085993269104826,
        ),
    ];
    let mut simplex = VoronoiSimplex::new();
    simplex.reset(CsoPoint::single_point(vertices[0]));
    for vertex in &vertices[1..] {
        assert!(simplex.add_point(CsoPoint::single_point(*vertex)));
    }
    let projection = simplex.project_origin_and_reduce();
    assert!(simplex.dimension() < 3);
    assert!(
        projection.length() > 2e-7 && projection.length() < 1e-6,
        "projection={projection:?} simplex={simplex:?}"
    );
    let mut witness = Vector::ZERO;
    for i in 0..=simplex.dimension() {
        witness += simplex.point(i).point * simplex.proj_coord(i);
    }
    assert!((projection - witness).length() < 1e-12);
}
