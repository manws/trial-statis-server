import asyncio
from typing import List
from app.core.mysql_pool import MySQLPool
from app.core.redis_pool import RedisPool
from app.utils.LoggerHelper import LoggerHelper

class AppInitializer:
    """应用初始化管理器"""
    
    def __init__(self):
        self.initialized_components: List[str] = []
        self.is_initialized = False
    
    async def initialize_all(self):
        """初始化所有应用组件"""
        if self.is_initialized:
            LoggerHelper.warning("应用组件已初始化，跳过重复初始化")
            return
            
        LoggerHelper.info("开始初始化应用组件...")
        
        try:
            # 1. 初始化Redis连接池
            LoggerHelper.info("1. 初始化Redis连接池...")
            await RedisPool.init_pool()
            self.initialized_components.append("redis")
            LoggerHelper.info("Redis连接池初始化完成")
            
            # 2. 初始化MySQL连接池（会在首次使用时自动初始化）
            LoggerHelper.info("2. MySQL连接池将在首次使用时自动初始化")
            self.initialized_components.append("mysql")
            
            # 3. 其他初始化工作可以在这里添加
            # 例如：加载配置、初始化缓存、预热数据等
            
            self.is_initialized = True
            LoggerHelper.info("所有应用组件初始化完成")
            
        except Exception as e:
            LoggerHelper.error(f"应用初始化失败: {e}")
            # 初始化失败时清理已初始化的组件
            await self.cleanup()
            raise
    
    async def cleanup(self):
        """清理所有已初始化的组件"""
        if not self.initialized_components:
            LoggerHelper.info("没有需要清理的组件")
            return
            
        LoggerHelper.info("开始清理应用组件...")
        cleanup_start = asyncio.get_event_loop().time()
        
        # 按相反顺序清理组件
        for component in reversed(self.initialized_components):
            try:
                if component == "mysql":
                    LoggerHelper.info("关闭MySQL连接池...")
                    await MySQLPool.close_pool()
                    LoggerHelper.info("MySQL连接池已关闭")
                    
                elif component == "redis":
                    LoggerHelper.info("关闭Redis连接池...")
                    await RedisPool.close_pool()
                    LoggerHelper.info("Redis连接池已关闭")
                    
            except Exception as e:
                LoggerHelper.error(f"关闭 {component} 组件时出错: {e}")
        
        self.initialized_components.clear()
        self.is_initialized = False
        
        cleanup_duration = asyncio.get_event_loop().time() - cleanup_start
        LoggerHelper.info(f"组件清理完成，耗时: {cleanup_duration:.2f}秒")

# 全局初始化器实例
app_initializer = AppInitializer()

# 便捷函数
async def init_app_params():
    """初始化应用参数的便捷函数"""
    await app_initializer.initialize_all()

async def cleanup_app_params():
    """清理应用参数的便捷函数"""
    await app_initializer.cleanup()

__all__ = ["AppInitializer", "app_initializer", "init_app_params", "cleanup_app_params"]