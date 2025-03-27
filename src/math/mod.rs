use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// A wrapper type for Vector3 that implements Serialize and Deserialize
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn to_vector3(&self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }

    pub fn from_vector3(v: &Vector3<f64>) -> Self {
        Self {
            x: v[0],
            y: v[1],
            z: v[2],
        }
    }
}

impl From<Vector3<f64>> for Vec3 {
    fn from(v: Vector3<f64>) -> Self {
        Self::from_vector3(&v)
    }
}

impl From<Vec3> for Vector3<f64> {
    fn from(v: Vec3) -> Self {
        v.to_vector3()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json;

    #[test]
    fn test_vec3_serialization() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let json = serde_json::to_string(&v).unwrap();
        let v2: Vec3 = serde_json::from_str(&json).unwrap();
        assert_eq!(v.x, v2.x);
        assert_eq!(v.y, v2.y);
        assert_eq!(v.z, v2.z);
    }

    #[test]
    fn test_vec3_conversion() {
        let v1 = Vector3::new(1.0, 2.0, 3.0);
        let v2 = Vec3::from_vector3(&v1);
        let v3 = v2.to_vector3();
        assert_eq!(v1, v3);
    }
}
