# Cellular Coverage Mapper

Matlab scripts for investigating system-level coverage for cellular networks, possibly with data relay drones, over a large geographic area. More information can be found in our paper:

>**Zhang, Y.**, Arakawa, T., Krogmeier, J.V., Anderson, C.R., Love, D.J. and Buckmaster, D.R., 2020, June. **Large-scale cellular coverage analyses for UAV data relay via channel modeling**. In *2020 IEEE International Conference on Communications (ICC)* (pp. 1-6). IEEE. DOI: [10.1109/ICC40277.2020.9149403](https://doi.org/10.1109/ICC40277.2020.9149403). [**\[Virtual presentation\]**](https://yaguangzhang.github.io/files/ICC2020_WC17_CellCoverageSimulationForDrones.mp4)

## Notes

For large simulations which take more than a couple of days to finish, Matlab `parfor` may freeze: the program seems to be running but the CPU usage is extremely low. If this happens, one can simply stop and restart the simulation, since the program is able to resume from progress log files.

## Contact

* **Yaguang Zhang** | *Purdue University* | Email: ygzhang@purdue.edu

## License

This project is licensed under the MIT License.
