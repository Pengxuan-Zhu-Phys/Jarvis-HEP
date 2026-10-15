# Jarvis-HEP 2.0.15 发行候选

1. **Roadmap phase**：2.0.15 可靠性修复及完整发行候选交付。
2. **Affected modules**：cardload、Worker/Factory/supervisor、monitor、CLI、nested sampler、发行工具；Jarvis-Operas help。
3. **Public interfaces**：命令及任务 YAML 保持兼容；修正表达式诊断、失败退出、monitor JSON 和缺失数据表达。发行工具新增 `--operas-dist`。
4. **Main risks**：表达式编译无法预测运行时数值域错误；Worker 启动默认等待上限 120 秒；macOS/Windows CPU 是活进程估算；nested 未报告终止准则时不推测收敛。
5. **Verification performed**：定向回归、spawn/Redis/恢复/吞吐验收、独立源码 fast gate；最终完整 wheel/sdist 干净安装验证通过；源码 fast gate 1272 passed，安装包回归 22 passed。

状态：**发行候选就绪，未发布**。最终发布由用户执行。

## 修复

- `nosuch.f(z)` 等不存在的 Operas 函数在 validation 阶段返回 `JV2-EXP-001`、具体 YAML 路径和纠正提示，不调用数值物理函数。Mapper 保持已有诊断。
- Worker 初始化错误通过 spawn pipe 反馈到控制进程；正常扫描前有界等待 ready/error，失败时停止进程并执行既有清理。运行期崩溃恢复与样本失败仍保留原语义。
- `Jarvis monitor --once --json R1` / 按名称选择输出一个 JSON 对象；无 ref 保留列表。读取 OS Worker 清单和对应 Redis boards，报告实际采集时间与 PID，未知数量为 null。连接/选择错误同样返回 JSON 与非零状态。
- task card 交给 plot 时明确要求带 `DataSet` 的绘图 YAML，提示 `Jarvis gen-plot-yaml TASK.yaml`。
- HEP intent flag 放在子命令前后都能报告 `conflicting CLI intents`；依赖 CLI 参数仍完整转发。
- Jarvis-Operas 1.4.2 用 Rich `Styled` 保留帮助 renderable，修复 `--force` 对象字符串。
- macOS/Windows scan CPU 按 PID+创建时间累计实际观察到的增量，避免退出进程导致旧读数冻结；标记为活进程估算，可能漏掉两次刷新之间结束的进程。
- nested 正常返回只报告 `engine returned`；未获得引擎明确原因时不把预算或自定义停止说成收敛。
- 吞吐验收预热完整管线并测量 16 个点，保持原 1.25 倍速度要求。

## 发行范围

包含本分支已跟踪的 Sampler plug-in/项目打包功能及本次修复，依赖 `Jarvis-HEP-Sampler>=0.1.0,<0.2` 和 `Jarvis-Operas>=1.4.2`。完整依赖一次安装，没有额外的精简包或 editable 安装。

未验收 Diver 源码、测试、两项 package-data 改动及其 F1/F2/F5 问题保留在工作区，不进入此次候选。Operas 的无关 hts 文档修改也未纳入依赖候选。工作区外部 Examples 卡片的非 ASCII 描述不属于 HEP 发行内容；独立候选测试不依赖该未发布卡片。

## 用户最终发布步骤

1. 审阅本地修复和台账；仅提交已验收文件，保留未完成 Diver 和无关 hts 改动。所有新源码、测试和维护文档必须进入 Git，`jhrel` 只快照 tracked files。
2. 先发布 Sampler SDK **0.1.0**，再发布 Jarvis-Operas **1.4.2**。两者候选文件位于 HEP 的 `dist/2.0.15/dependencies/`；其 SHA256 与验证结果随候选交付。
3. 依赖公开后，从已验收 release commit 的干净 checkout 使用新的空目录再执行正常发布前验证：

   ```bash
   ./jhrel --output dist/2.0.15-public-check
   ```

   本地依赖候选验收可用 `--sampler-dist dist/2.0.15/dependencies/sampler --operas-dist dist/2.0.15/dependencies/operas`。这两个参数仅替换依赖来源，仍安装每项声明依赖。
4. 用户将已验收 release commit 合入 `master`，自行创建 `v2.0.15` GitHub Release；现有 Trusted Publishing workflow 构建、验证并发布到 PyPI。正式发布后核验公开 wheel 版本与 `Jarvis -v`。本轮代理不创建 tag、推送或上传。

版本是否已经存在，发布时应再次检查。候选就绪与公开发行是两个独立状态。
