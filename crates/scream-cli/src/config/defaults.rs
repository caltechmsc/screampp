pub struct DefaultsConfig {
    pub s_factor: f64,
    pub forcefield: String,
    pub delta_params: String,
    pub topology_registry: String,
    pub rotamer_library: String,
    pub num_solutions: usize,
    pub max_iterations: usize,
    pub include_input_conformation: bool,
    pub final_refinement_iterations: usize,
    pub energy_threshold: f64,
    pub patience_iterations: usize,
}

impl Default for DefaultsConfig {
    fn default() -> Self {
        Self {
            s_factor: 1.1,
            forcefield: "exp-6@0.4".to_string(),
            delta_params: "rmsd-1.0".to_string(),
            topology_registry: "default".to_string(),
            rotamer_library: "charmm@rmsd-1.0".to_string(),
            include_input_conformation: true,
            num_solutions: 1,
            max_iterations: 100,
            final_refinement_iterations: 2,
            energy_threshold: 0.01,
            patience_iterations: 5,
        }
    }
}
