/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include <time.h>
#include <stdlib.h>
#include "etCore.h"
#include "etExample.h"
#include "myOptparser.h"


INT usage(const CHAR *program)
{
    INT a,c;
    print_err("usage: %s [options]\n", program);
    print_err("\t-h, --help\tprint usage and exit\n");
    print_err("\t-s SOURCE, --source=SOURCE\n\t\ttype of source terms [0]\n");
    print_err("\t\t0: sourceless\n");

    // print_err("\t-m MASS, --mass=MASS\n\t\tmass of central Kerr BH [%g]\n", DEFAULT_MASS);
    print_err("\t-a SPIN, --spin=SPIN\n\t\tspin of central Kerr BH [%g]\n", DEFAULT_SPIN);

    print_err("\t-r MASS-RATIO, --mass-ratio=MASS-RATIO\n\t\tmass ratio between secondary mass and central Kerr BH [%g]\n", DEFAULT_ETA);

    print_err("\t-R STEP-RATIO, --step-ratio=STEP-RATIO\n\t\tratio between time step and radius step [%g]\n", DEFAULT_STEP_RATIO);
    print_err("\t-T TIME-MAX, --time-max=TIME-MAX\n\t\tmax evolution time [%gM]\n", DEFAULT_TIME_MAX);
    print_err("\t-S SHOT-STEP, --shot-step=TIME-MAX\n\t\ttime interval of save [%gM]\n", 50);

    print_err("\t-p PREFIX, --prefix=PREFIX\n\t\toutput path [TeukDefault]\n");
    print_err("\t-v VERSION, --version=VERSION\tVersion of code[1]\n\t\t\n");
    print_err("\t-E MODE --mode=MODE        \n\t\tThe output mode.\n");
    print_err("\t-l LOG_LEVEL, --log-level=LOG_LEVEL\n\t\tDebug level for logging [0]\n");
    print_err("\t-o LOG_FILE, --log-file=LOG_FILE\n\t\tWill dump loggings to here.\n");
    print_err("\n\n");
    return CEV_SUCCESS;
}

CoreParams parseargs(INT argc, CHAR **argv, CtrlParams *ctpms)
{
    CoreParams p;
    p.mass = DEFAULT_MASS;
    p.spin = DEFAULT_SPIN;
    p.eta = DEFAULT_ETA;
    p.mode = 0;
    p.step_ratio = DEFAULT_STEP_RATIO;
    p.t_max = DEFAULT_TIME_MAX;
    p.source_type = 0;
    p.rt0 = 20.;
    p.theta0 = CST_PI/4.;
    p.mode = 0;
    p.shot_step = 50;
    SET_CODE_VERSION(1);
    extern CHAR *EXT_optarg;
    extern INT EXT_optind;
    ctpms->level = 1;
    strncpy(ctpms->flog, "None", STR_COMM_SIZE);
    strncpy(p.prefix, "TeukDefault", STR_COMM_SIZE);
    OPTION long_options[] = {
        {"help", opt_no_argument, 0, 'h'},

        {"source", opt_required_argument, 0, 's'},
        {"spin", opt_required_argument, 0, 'a'},
        {"mass-ratio", opt_required_argument, 0, 'r'},

        {"step-ratio", opt_required_argument, 0, 'R'},
        {"time-max", opt_required_argument, 0, 'T'},
        {"shot-step", opt_required_argument, 0, 'S'},

        {"prefix", opt_required_argument , 0, 'p'},

        {"mode", opt_required_argument, 0, 'E'},
        {"version", opt_required_argument, 0, 'v'},
        {"log-level", opt_required_argument , 0, 'l'},
        {"log-file", opt_required_argument , 0, 'o'},
        {0, 0, 0, 0}
    };

    CHAR args[] =
    "h:s:m:a:r:R:T:p:E:v:l:o";
    while (1)
    {
        INT option_index = 0;
        INT c;
        c = getopt_long_only(argc, argv, args, long_options, &option_index);
        if (c == -1)
            break;
        switch (c)
        {
            case 0:
                if (long_options[option_index].flag)
                    break;
                else
                {
                    print_err("error parsing option %s with argument %s\n", long_options[option_index].name, EXT_optarg);
                    exit(1);
                }
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'v':
                SET_CODE_VERSION(atoi(EXT_optarg));
                break;
            case 's':
                p.source_type = atoi(EXT_optarg);
                break;
            case 'a':
                p.spin = atof(EXT_optarg);
                break;
            case 'r':
                p.eta = atof(EXT_optarg);
                break;
            case 'R':
                p.step_ratio = atof(EXT_optarg);
                break;
            case 'T':
                p.t_max = atof(EXT_optarg);
                break;
            case 'S':
                p.shot_step = atoi(EXT_optarg);
                break;
            case 'p':
                strncpy(p.prefix, EXT_optarg, STR_COMM_SIZE);
                break;
            case 'E':
                p.mode = atoi(EXT_optarg);
                break;
            case 'l':
                ctpms->level = atoi(EXT_optarg);
                break;
            case 'o':
                strncpy(ctpms->flog, EXT_optarg, STR_COMM_SIZE);
                break;
            default:
                print_err("unknown error while parsing options\n");
                exit(1);
        }
    }
    if (EXT_optind < argc)
    {
        print_err("extraneous command line arguments:\n");
        while (EXT_optind < argc)
            print_err("%s\n", argv[EXT_optind++]);
        exit(1);
    }
    return p;
}


int main(int  argc, char **argv)
{
    CoreParams p;
    INT status;
    CtrlParams ctpms;
    p = parseargs(argc, argv, &ctpms);
    if (strcmp(ctpms.flog, "None")!=0)
    {
        LOG_SetPrintLogPlaceFlag(1);
    }
    LOG_SetPrintDebugLogFlag(ctpms.level);
    LOG_Init(ctpms.flog, 40960);
    PRINT_LOG_INFO(LOG_INFO, "CodeVersionFlag = %d\n", CODE_VERSION);
    PRINT_LOG_INFO(LOG_INFO, "Output PATH: %s", p.prefix);
    status = cmd_mkdir(p.prefix);
    // status = Example_String_Vibration(&p);
    status = Solve_Teukolsky_TD(&p);
    if (status != CEV_SUCCESS)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "THE PROGRAM FAILED.");
        goto EXIT;
    }
    PRINT_LOG_INFO(LOG_INFO, "PROGRAM END.");
EXIT:
    CheckMemoryLeak();
    return 0;
}
