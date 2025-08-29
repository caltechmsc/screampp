use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use super::context::OptimizationContext;
use super::error::EngineError;
use super::placement;
use std::collections::HashMap;

pub struct SystemView<'a, 'ctx, C>
where
    C: super::context::ProvidesResidueSelections + Sync,
{
    pub system: &'a mut MolecularSystem,
    pub context: &'ctx OptimizationContext<'ctx, C>,
    pub current_rotamers: &'a mut HashMap<ResidueId, usize>,
}

impl<'a, 'ctx, C> SystemView<'a, 'ctx, C>
where
    C: super::context::ProvidesResidueSelections + Sync,
{
    pub fn new(
        system: &'a mut MolecularSystem,
        context: &'ctx OptimizationContext<'ctx, C>,
        current_rotamers: &'a mut HashMap<ResidueId, usize>,
    ) -> Self {
        Self {
            system,
            context,
            current_rotamers,
        }
    }

    pub fn apply_move(
        &mut self,
        res_id: ResidueId,
        new_rotamer_idx: usize,
    ) -> Result<(), EngineError> {
        self.place_rotamer(res_id, new_rotamer_idx)?;
        self.current_rotamers.insert(res_id, new_rotamer_idx);
        Ok(())
    }

    fn place_rotamer(&mut self, res_id: ResidueId, rotamer_idx: usize) -> Result<(), EngineError> {
        let residue = self.system.residue(res_id).unwrap();
        let res_type = residue.residue_type.unwrap();
        let rotamer = &self
            .context
            .rotamer_library
            .get_rotamers_for(res_type)
            .unwrap()[rotamer_idx];
        let res_name = res_type.to_three_letter();
        let topology = self.context.topology_registry.get(res_name).unwrap();

        placement::place_rotamer_on_system(self.system, res_id, rotamer, topology)
    }

    pub fn transaction<F, R>(
        &mut self,
        res_id_to_modify: ResidueId,
        action: F,
    ) -> Result<R, EngineError>
    where
        F: FnOnce(&mut Self) -> Result<R, EngineError>,
    {
        // 1. Record the original state for the residue.
        let original_rotamer_idx =
            *self
                .current_rotamers
                .get(&res_id_to_modify)
                .ok_or_else(|| {
                    EngineError::Internal(format!(
                        "Residue {:?} not found in rotamer tracking map",
                        res_id_to_modify
                    ))
                })?;

        // 2. Execute the action.
        let result = action(self)?;

        // 3. Check if the residue was actually modified during the action.
        let current_rotamer_idx = *self.current_rotamers.get(&res_id_to_modify).unwrap();

        // 4. If it was modified, revert it back to the original state.
        if current_rotamer_idx != original_rotamer_idx {
            self.place_rotamer(res_id_to_modify, original_rotamer_idx)?;
            self.current_rotamers
                .insert(res_id_to_modify, original_rotamer_idx);
        }

        Ok(result)
    }
    pub fn transaction_doublet<F, R>(
        &mut self,
        res_a_id: ResidueId,
        res_b_id: ResidueId,
        action: F,
    ) -> Result<R, EngineError>
    where
        F: FnOnce(&mut Self) -> Result<R, EngineError>,
    {
        // 1. Record the original states for both residues.
        let original_rotamer_idx_a = *self.current_rotamers.get(&res_a_id).ok_or_else(|| {
            EngineError::Internal(format!(
                "Residue {:?} not found in rotamer tracking map",
                res_a_id
            ))
        })?;

        let original_rotamer_idx_b = *self.current_rotamers.get(&res_b_id).ok_or_else(|| {
            EngineError::Internal(format!(
                "Residue {:?} not found in rotamer tracking map",
                res_b_id
            ))
        })?;

        // 2. Execute the action.
        let result = action(self)?;

        // 3. Check if the residues were actually modified during the action.
        let current_rotamer_idx_a = *self.current_rotamers.get(&res_a_id).unwrap();
        let current_rotamer_idx_b = *self.current_rotamers.get(&res_b_id).unwrap();

        // 4. If they were modified, revert them back to the original states.
        if current_rotamer_idx_a != original_rotamer_idx_a {
            self.place_rotamer(res_a_id, original_rotamer_idx_a)?;
            self.current_rotamers
                .insert(res_a_id, original_rotamer_idx_a);
        }

        if current_rotamer_idx_b != original_rotamer_idx_b {
            self.place_rotamer(res_b_id, original_rotamer_idx_b)?;
            self.current_rotamers
                .insert(res_b_id, original_rotamer_idx_b);
        }

        Ok(result)
    }
}
