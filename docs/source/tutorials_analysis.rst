Analysis and Error Estimation
=============================

Analysis on partitions
----------------------

To analyse the statistics of a ``HEALPix`` map -- ``map1``, for each partition simply isolate the pixels
of each partition and then run the statistics on those map points.

.. code-block:: python

    for partID in range(1, Npartition+1):

        # for map based partitions
        cond = np.where(partitionmap == partID)[0]

        map_point_in_partition = map1[cond]

        """
        INSERT ANALYSIS ON MAP POINTS HERE 
        """

        # or create a new map with only the partitioned points on the map

        map_partition = np.zeros(hp.nside2npix(nside))
        map_partition[cond] = map1[cond]

        """
        INSERT ANALYSIS ON PARTITION MAP
        """

To analyse the statistics of the statistics ``stats`` of a set of points we simply do the following

.. code-block:: python

    for partID in range(1, Npartition+1):

        # for map based partitions
        cond = np.where(partitionID == partID)[0]

        point_in_partition = stats[cond]

        """
        INSERT ANALYSIS ON POINTS HERE 
        """

Jackknife error estimation
--------------------------

For estimating jackknife based errors we will want all the points in a map that are not in a 
chosen partition, this can be done by running:

.. code-block:: python

    stats_jk = []

    for partID in range(1, Npartition+1):

        # for map based partitions
        cond = np.where((partitionmap != partID) & (partitionmap != 0))[0]

        map_point_in_partition = map1[cond]

        """
        INSERT ANALYSIS ON MAP POINTS HERE 
        to create stats_on_partition data vector
        """

        # or create a new map with only the partitioned points on the map

        map_partition = np.zeros(hp.nside2npix(nside))
        map_partition[cond] = map1[cond]

        """
        INSERT ANALYSIS ON PARTITION MAP 
        to create stats_on_partition data vector
        """

        # store the statistics data vector
        stats_jk.append(stats_on_partition)
    
    stats_jk = np.array(stats_jk)

    # jackknife mean and errors can be computed by
    stats_mean = np.mean(stats_jk, axis=0)

    jk_prefactor = (Npartition - 1)*(Npartition - 1)/Npartition
    stats_std = np.sqrt(jk_prefactor)*np.std(stats_jk, axis=0)

    # while the covariance can be computed by
    stats_cov = jk_prefactor * np.cov(stats_jk.T)


.. note:: 

    We must be careful to also remove the masked region of the map which is given a ``partitionID = 0``.

For points with a statistiss ``stats`` we do something very similar except we do not need to worry about 
the masked region here.

.. code-block:: python

    stats_jk = []

    for partID in range(1, Npartition+1):

        # for map based partitions
        cond = np.where(partitionID != partID)[0]

        point_in_partition = stats[cond]

        """
        INSERT ANALYSIS ON POINTS HERE 
        to create stats_on_partition data vector
        """

        # store the statistics data vector
        stats_jk.append(stats_on_partition)
    
    stats_jk = np.array(stats_jk)

    # jackknife mean and errors can be computed by
    stats_mean = np.mean(stats_jk, axis=0)
    
    jk_prefactor = (Npartition - 1)*(Npartition - 1)/Npartition
    stats_std = np.sqrt(jk_prefactor)*np.std(stats_jk, axis=0)

    # while the covariance can be computed by
    stats_cov = jk_prefactor * np.cov(stats_jk.T)