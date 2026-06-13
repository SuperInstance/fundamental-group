# Fundamental Group Calculator

**A Rust library for computing fundamental groups π₁(X, x₀) of topological spaces**, implementing path composition, group presentations, the Seifert-van Kampen theorem, and winding number computation.

## Why It Matters

The fundamental group is the most basic algebraic invariant in algebraic topology — it classifies loops in a space up to continuous deformation. It distinguishes a circle (π₁ = ℤ) from a sphere (π₁ = 0), detects holes in a surface, and underlies covering space theory. In robotics, fundamental groups appear in motion planning (homotopy classes of paths). In physics, they classify topological defects in condensed matter (vortices, monopoles). This library provides both the algebraic machinery (group presentations with generators and relations) and geometric tools (winding numbers, path composition).

## How It Works

The library models two layers of topology: **geometric** and **algebraic**.

**Geometric layer**: `Path` is a sequence of points in ℝⁿ. Paths can be composed end-to-end (if the endpoint of one matches the start of another), reversed, and checked for being a loop. The `Homotopy` module computes winding numbers — the net number of times a closed curve encircles the origin — by accumulating angle changes using `atan2` and normalizing to `[-π, π]`.

**Algebraic layer**: `GroupPresentation` represents a group as `⟨generators | relations⟩`. Generators are named symbols with optional inverses. Relations are words (sequences of generators). The library provides constructors for common spaces: `of_sphere(n)` returns the trivial group for n ≥ 2 and ℤ for n = 1; `of_torus(n)` returns the free abelian group ℤⁿ; `of_rp2()` returns ℤ/2ℤ; `of_klein_bottle()` returns `⟨a,b | aba⁻¹b⟩`. The `van_kampen` method applies the Seifert-van Kampen theorem to compute π₁ of a union from the fundamental groups of its parts. The `abelianize` method converts any presentation to its abelianization (the first homology group H₁) by tracking generator cancellation in commutator relations.

## Quick Start

```rust
use fundamental_group::*;

fn main() {
    // Fundamental group of common spaces
    let s1 = FundamentalGroup::of_sphere(1);
    println!("π₁(S¹) = {} generators", s1.generators.len()); // 1 (ℤ)

    let s2 = FundamentalGroup::of_sphere(2);
    println!("π₁(S²) = {} generators", s2.generators.len()); // 0 (trivial)

    let torus = FundamentalGroup::of_torus(2);
    println!("π₁(T²) = {} generators", torus.generators.len()); // 2 (ℤ²)

    // Klein bottle: ⟨a,b | aba⁻¹b⟩
    let klein = FundamentalGroup::of_klein_bottle();
    println!("π₁(Klein) = {{{}}}", klein.generators.join(", "));

    // Compute winding number of a loop around the origin
    let loop_path = Path::new(vec![
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![-1.0, 0.0],
        vec![0.0, -1.0],
        vec![1.0, 0.0],
    ]);
    println!("Winding number: {}", Homotopy::winding_number(&loop_path)); // 1
}
```

## API

| Type / Function | Description |
|---|---|
| `Path::new(points)` | Construct a path from a sequence of coordinates |
| `Path::compose(other)` | Concatenate two compatible paths |
| `Path::is_loop()` | Check if path starts and ends at the same point |
| `GroupPresentation` | `⟨generators | relations⟩` presentation |
| `FundamentalGroup::of_sphere(n)` | π₁(Sⁿ) — ℤ for n=1, trivial for n≥2 |
| `FundamentalGroup::of_torus(n)` | π₁(Tⁿ) = ℤⁿ |
| `FundamentalGroup::of_rp2()` | π₁(RP²) = ℤ/2ℤ |
| `FundamentalGroup::of_klein_bottle()` | π₁(Klein) = `⟨a,b | aba⁻¹b⟩` |
| `van_kampen(π_U, π_V, π_U∩V)` | Seifert-van Kampen theorem |
| `Homotopy::winding_number(path)` | Winding number of a planar loop |

## Architecture Notes

Part of the SuperInstance topology suite. Companion crates: `floer-homology`, `hodge-theory`, `handle-decomposition`, `hurewicz-homomorphism`. See the [Architecture Guide](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## License

MIT
