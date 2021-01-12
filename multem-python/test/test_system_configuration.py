import multem
import pickle


def test_system_configuration():

    system_conf = multem.System_Configuration()
    system_conf.device = "device"
    system_conf.precision = "float"
    system_conf.cpu_ncores = 1
    system_conf.cpu_nthread = 2
    system_conf.gpu_device = 3
    system_conf.gpu_nstream = 4

    def check():
        assert system_conf.device.name == "device"
        assert system_conf.precision.name == "float"
        assert system_conf.cpu_ncores == 1
        assert system_conf.cpu_nthread == 2
        assert system_conf.gpu_device == 3
        assert system_conf.gpu_nstream == 4

    check()

    system_conf = pickle.loads(pickle.dumps(system_conf))

    check()
