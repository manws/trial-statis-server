from fastapi import APIRouter
from .base_controller import router as base_router
from .anova_controller import router as anova_router
from .chi_controller import router as chi_router
from .clu_controller import router as clu_router
from .cor_controller import router as cor_router
from .cq_controller import router as cq_router
from .hv_controller import router as hv_router
from .reg_controller import router as reg_router
from .rs_controller import router as rs_router
from .runs_controller import router as runs_router
from .sur_controller import router as sur_router
from .t_controller import router as t_router
from .z_controller import router as z_router

router = APIRouter(prefix="/stats")

router.include_router(base_router)
router.include_router(anova_router)
router.include_router(chi_router)
router.include_router(clu_router)
router.include_router(cor_router)
router.include_router(cq_router)
router.include_router(hv_router)
router.include_router(reg_router)
router.include_router(rs_router)
router.include_router(runs_router)
router.include_router(sur_router)
router.include_router(t_router)
router.include_router(z_router)

# expose for parent imports
stats_router = router