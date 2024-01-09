def betti_number(self, p: int, group=None) -> int:
        """
        Gets the betti numbers of the simplicial complex for the given dimension p.

        Args:
        p (int): dimension
        group: group

        Returns:
        int: betti_number
        """
        if group is None:
            mp = self.smith_normal_form_z(self.boundary_matrix(p, group))
            mp_1 = self.smith_normal_form_z(self.boundary_matrix(p + 1, group))
        else:
            mp = self.smith_normal_form(self.boundary_matrix(p, group), group=group)
            mp_1 = self.smith_normal_form(self.boundary_matrix(p + 1, group), group=group)

        # Number of columns of zeros
        dim_zp = len([_ for x in np.transpose(mp) if sum(x) == 0])

        # Number of rows with ones
        dim_bp = len([_ for x in mp_1 if sum(x) != 0])

        return dim_zp - dim_bp