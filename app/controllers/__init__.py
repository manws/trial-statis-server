from fastapi import APIRouter

from .stats import stats_router

api_router = APIRouter(prefix="/api/v1")

# 包含所有子路由
api_router.include_router(stats_router)
