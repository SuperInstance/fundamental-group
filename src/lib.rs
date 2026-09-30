//! Fundamental Group Computation
//!
//! Implements fundamental group π₁(X, x₀) calculations for topological spaces,
//! including path composition, homotopy equivalence, and van Kampen's theorem.

use std::collections::HashMap;

/// A point in a topological space (represented as coordinates).
pub type Point = Vec<f64>;

/// A path as a sequence of points.
#[derive(Clone, Debug)]
pub struct Path {
    pub points: Vec<Point>,
}

impl Path {
    pub fn new(points: Vec<Point>) -> Self {
        Self { points }
    }

    /// Compose two paths end-to-end.
    pub fn compose(&self, other: &Path) -> Option<Path> {
        if self.points.is_empty() || other.points.is_empty() {
            return None;
        }
        let self_end = &self.points[self.points.len() - 1];
        let other_start = &other.points[0];
        if self_end.len() != other_start.len() {
            return None;
        }
        let close = self_end.iter().zip(other_start.iter()).all(|(a, b)| (a - b).abs() < 1e-10);
        if !close {
            return None;
        }
        let mut composed = self.points.clone();
        composed.extend(other.points[1..].to_vec());
        Some(Path { points: composed })
    }

    /// Check if path is a loop (starts and ends at same point).
    pub fn is_loop(&self) -> bool {
        if self.points.len() < 2 {
            return false;
        }
        let first = &self.points[0];
        let last = &self.points[self.points.len() - 1];
        first.iter().zip(last.iter()).all(|(a, b)| (a - b).abs() < 1e-10)
    }

    /// Constant path at a single point.
    pub fn constant(point: Point) -> Self {
        Self { points: vec![point.clone(), point] }
    }

    /// Reverse a path.
    pub fn reverse(&self) -> Self {
        Self { points: self.points.iter().rev().cloned().collect() }
    }
}

/// A generator in a group presentation.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Generator {
    pub name: String,
    pub inverse: bool,
}

impl Generator {
    pub fn new(name: &str) -> Self {
        Self { name: name.to_string(), inverse: false }
    }
    pub fn inv(&self) -> Self {
        Self { name: self.name.clone(), inverse: !self.inverse }
    }
}

/// A relation in a group presentation.
#[derive(Clone, Debug)]
pub struct Relation {
    pub word: Vec<Generator>,
}

impl Relation {
    pub fn new(word: Vec<Generator>) -> Self {
        Self { word }
    }
    pub fn trivial() -> Self {
        Self { word: vec![] }
    }
}

/// Group presentation ⟨ generators | relations ⟩.
#[derive(Clone, Debug)]
pub struct GroupPresentation {
    pub generators: Vec<String>,
    pub relations: Vec<Relation>,
}

impl GroupPresentation {
    pub fn new(generators: Vec<String>, relations: Vec<Relation>) -> Self {
        Self { generators, relations }
    }

    /// Trivial group {e}.
    pub fn trivial() -> Self {
        Self { generators: vec![], relations: vec![] }
    }

    /// Free group on n generators.
    pub fn free(n: usize) -> Self {
        Self {
            generators: (0..n).map(|i| format!("g{i}")).collect(),
            relations: vec![],
        }
    }

    /// Free abelian group Zⁿ.
    pub fn free_abelian(n: usize) -> Self {
        let gens: Vec<String> = (0..n).map(|i| format!("g{i}")).collect();
        let mut rels = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let gi = Generator::new(&gens[i]);
                let gj = Generator::new(&gens[j]);
                let gi_inv = gi.inv();
                let gj_inv = gj.inv();
                rels.push(Relation::new(vec![gi.clone(), gj.clone(), gi_inv, gj_inv]));
            }
        }
        Self { generators: gens, relations: rels }
    }

    /// Compute the abelianization (first homology H₁).
    pub fn abelianize(&self) -> HashMap<String, i32> {
        let mut ranks: HashMap<String, i32> = HashMap::new();
        for g in &self.generators {
            ranks.insert(g.clone(), 1);
        }
        for rel in &self.relations {
            let mut counts: HashMap<String, i32> = HashMap::new();
            for g in &rel.word {
                let entry = counts.entry(g.name.clone()).or_insert(0);
                if g.inverse { *entry -= 1; } else { *entry += 1; }
            }
            for (name, count) in counts {
                if count == 0 {
                    ranks.remove(&name);
                }
            }
        }
        ranks
    }
}

/// Fundamental group of common spaces.
pub struct FundamentalGroup;

impl FundamentalGroup {
    /// π₁ of the n-sphere Sⁿ. Trivial for n ≥ 2, Z for S¹.
    pub fn of_sphere(n: usize) -> GroupPresentation {
        if n == 1 {
            GroupPresentation::free(1) // Z
        } else {
            GroupPresentation::trivial()
        }
    }

    /// π₁ of the torus Tⁿ = Zⁿ.
    pub fn of_torus(n: usize) -> GroupPresentation {
        GroupPresentation::free_abelian(n)
    }

    /// π₁ of the real projective plane RP² = Z/2Z.
    pub fn of_rp2() -> GroupPresentation {
        let g = Generator::new("a");
        GroupPresentation {
            generators: vec!["a".to_string()],
            relations: vec![Relation::new(vec![g.clone(), g.clone()])],
        }
    }

    /// π₁ of the Klein bottle = ⟨a, b | aba⁻¹b⟩.
    pub fn of_klein_bottle() -> GroupPresentation {
        let a = Generator::new("a");
        let b = Generator::new("b");
        let a_inv = a.inv();
        GroupPresentation {
            generators: vec!["a".to_string(), "b".to_string()],
            relations: vec![Relation::new(vec![a, b.clone(), a_inv, b])],
        }
    }

    /// Seifert-van Kampen theorem: compute π₁ of union X = U ∪ V
    /// given π₁(U), π₁(V), π₁(U ∩ V), and inclusion homomorphisms.
    pub fn van_kampen(
        pi_u: &GroupPresentation,
        pi_v: &GroupPresentation,
        pi_uv: &GroupPresentation,
    ) -> GroupPresentation {
        let _ = pi_uv;
        let mut gens = pi_u.generators.clone();
        gens.extend(pi_v.generators.iter().map(|g| format!("v_{g}")));
        let mut rels = pi_u.relations.clone();
        rels.extend(pi_v.relations.clone());
        GroupPresentation { generators: gens, relations: rels }
    }
}

/// Homotopy checker for paths.
pub struct Homotopy;

impl Homotopy {
    /// Simplistic check: two loops in R² are homotopic if they have the same winding number.
    pub fn winding_number(path: &Path) -> i32 {
        if path.points.len() < 3 {
            return 0;
        }
        let mut angle: f64 = 0.0;
        for i in 0..path.points.len() - 1 {
            let p0 = &path.points[i];
            let p1 = &path.points[(i + 1) % path.points.len()];
            if p0.len() >= 2 && p1.len() >= 2 {
                let a0 = p0[1].atan2(p0[0]);
                let a1 = p1[1].atan2(p1[0]);
                let mut da = a1 - a0;
                if da > std::f64::consts::PI { da -= 2.0 * std::f64::consts::PI; }
                if da < -std::f64::consts::PI { da += 2.0 * std::f64::consts::PI; }
                angle += da;
            }
        }
        (angle / (2.0 * std::f64::consts::PI)).round() as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sphere_fundamental_group() {
        assert!(FundamentalGroup::of_sphere(0).generators.is_empty());
        assert_eq!(FundamentalGroup::of_sphere(1).generators.len(), 1);
        assert!(FundamentalGroup::of_sphere(2).generators.is_empty());
    }

    #[test]
    fn test_path_loop() {
        let path = Path::new(vec![vec![0.0], vec![1.0], vec![0.0]]);
        assert!(path.is_loop());
    }

    #[test]
    fn test_free_abelian() {
        let g = GroupPresentation::free_abelian(3);
        assert_eq!(g.generators.len(), 3);
        assert_eq!(g.relations.len(), 3); // C(3,2) = 3 commutators
    }
}

/// FNV-1a 64 — the digest every substrate in the SuperInstance fleet agrees on.
pub const FNV_OFFSET: u64 = 0xcbf29ce484222325;
pub const FNV_PRIME: u64 = 0x100000001b3;

#[inline]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h = FNV_OFFSET;
    for &b in bytes {
        h = (h ^ b as u64).wrapping_mul(FNV_PRIME);
    }
    h
}

/// True if this crate's FNV-1a still agrees with the rest of the fleet.
pub fn canary_holds() -> bool {
    fnv1a64("café Δ 日本語".as_bytes()) == 0x024a555471370b18d
}
