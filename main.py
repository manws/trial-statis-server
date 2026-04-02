import signal
import platform
import uvicorn
from contextlib import asynccontextmanager
from fastapi import FastAPI
from app.controllers import api_router
from app.utils.LoggerHelper import LoggerHelper
from app.core.app_init import app_initializer

@asynccontextmanager
async def program_life(app: FastAPI):
    """应用生命周期管理"""
    try:
        # 使用应用初始化管理器初始化所有组件
        await app_initializer.initialize_all()
        yield
        
    finally:
        # 关闭时清理资源
        LoggerHelper.info("应用关闭中，正在清理资源...")
        await app_initializer.cleanup()
        LoggerHelper.info("资源清理完成，应用已安全关闭")

# 创建FastAPI应用实例
app = FastAPI(
    title="Trial Statistics Server",
    description="临床试验统计服务API",
    version="1.0.0",
    lifespan=program_life 
)

# 注册路由
app.include_router(api_router)

def setup_signal_handlers():
    """设置跨平台信号处理器"""
    current_platform = platform.system().lower()
    
    def signal_handler(signum, frame):
        """通用信号处理器"""
        signal_name = signal.Signals(signum).name if hasattr(signal, 'Signals') else f"Signal {signum}"
        LoggerHelper.info(f"接收到 {signal_name} 信号 ({signum})，准备优雅关闭...")
    
    # 注册信号处理器
    if current_platform == "windows":
        # Windows环境
        try:
            signal.signal(signal.SIGINT, signal_handler)    # Ctrl+C
            signal.signal(signal.SIGBREAK, signal_handler)  # Ctrl+Break (Windows特有)
            LoggerHelper.info("Windows信号处理器注册完成 (SIGINT, SIGBREAK)")
        except AttributeError:
            # 某些Windows版本可能不支持SIGBREAK
            signal.signal(signal.SIGINT, signal_handler)
            LoggerHelper.info("Windows信号处理器注册完成 (SIGINT only)")
    else:
        # Unix/Linux环境 (包括Docker)
        try:
            signal.signal(signal.SIGTERM, signal_handler)   # DOCKER和Unix标准终止信号
            signal.signal(signal.SIGINT, signal_handler)    # Ctrl+C
            LoggerHelper.info("Unix信号处理器注册完成 (SIGTERM, SIGINT)")
        except Exception as e:
            LoggerHelper.warning(f"信号注册警告: {e}")
            # 至少注册SIGINT
            signal.signal(signal.SIGINT, signal_handler)
            LoggerHelper.info("基本信号处理器注册完成 (SIGINT)")

if __name__ == "__main__":
    LoggerHelper.info(f"Starting Trial Statistics Server on {platform.system()}...")
    
    # 设置跨平台信号处理器
    setup_signal_handlers()
    
    try:
        uvicorn.run(
            "main:app",
            host="0.0.0.0",
            port=8000,
            reload=True,  # 开发环境启用热重载
            log_level="info"
        )
    except KeyboardInterrupt:
        LoggerHelper.info("收到键盘中断，正在关闭服务...")
    except Exception as e:
        LoggerHelper.error(f"服务运行出错: {e}")
    finally:
        LoggerHelper.info("服务已停止")