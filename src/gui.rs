use gpui::{
    App, Bounds, ClipboardItem, Context, Entity, IntoElement, Render, SharedString, Subscription,
    Task, TitlebarOptions, Window, WindowBounds, WindowOptions, div, prelude::*, px, rgb, size,
};
use gpui_base::{Button, ColorTokens, Root, Theme, ThemeAppearance};
use gpui_graph::{Edge, GraphData, LabelMode, NetworkGraph, NetworkGraphEvent, Node, Selection};

use crate::explorer::{Field, Model, Snapshot};

const BG: u32 = 0x10151d;
const PANEL: u32 = 0x171e29;
const LINE: u32 = 0x2b3647;
const TEXT: u32 = 0xe4eaf2;
const MUTED: u32 = 0x95a4b9;
const ACCENT: u32 = 0x71dfbf;
const CHANGED: u32 = 0xffcd78;
const PAGE_SIZE: usize = 30;

pub fn run() {
    gpui_platform::application().run(|cx: &mut App| {
        gpui_base::init(cx);
        gpui_graph::init(cx);
        let theme = Theme::global_mut(cx);
        theme.appearance = ThemeAppearance::Dark;
        theme.tokens.colors = ColorTokens::dark();
        let bounds = Bounds::centered(None, size(px(1280.), px(800.)), cx);
        cx.open_window(
            WindowOptions {
                window_bounds: Some(WindowBounds::Windowed(bounds)),
                window_min_size: Some(size(px(1000.), px(640.))),
                titlebar: Some(TitlebarOptions {
                    title: Some("State graph explorer".into()),
                    ..Default::default()
                }),
                ..Default::default()
            },
            |window, cx| {
                let view = cx.new(Explorer::new);

                cx.new(|cx| Root::new(view, window, cx))
            },
        )
        .expect("Could not open the graph explorer window");
        cx.activate(true);
    });
}

struct Explorer {
    graph: Entity<NetworkGraph>,
    _subscription: Subscription,
    build_task: Option<Task<()>>,
    model: Option<Model>,
    snapshot: Option<Snapshot>,
    error: Option<String>,
    loading: bool,
    limit: usize,
    selected: usize,
    previous: Option<usize>,
    edge: Option<usize>,
    context_tab: bool,
    expanded_field: Option<String>,
    state_page: usize,
    neighborhood: bool,
    edge_labels: bool,
    fit_pending: bool,
}

impl Explorer {
    fn new(cx: &mut Context<Self>) -> Self {
        let graph = cx.new(|cx| {
            NetworkGraph::new(cx)
                .with_label_mode(LabelMode::Focus)
                .edge_labels(false)
        });
        let subscription = cx.subscribe(&graph, |this, _, event, cx| match event {
            NetworkGraphEvent::SelectionChanged(Selection::Node(node))
            | NetworkGraphEvent::NodeActivated(node) => {
                if let Ok(index) = node.id.parse::<usize>() {
                    this.select(index, false, cx);
                }
            }
            NetworkGraphEvent::SelectionChanged(Selection::Edge(edge)) => {
                this.edge = edge
                    .properties
                    .iter()
                    .find(|(k, _)| k.as_ref() == "index")
                    .and_then(|(_, v)| v.parse().ok());
                this.context_tab = false;
                this.expanded_field = None;
                cx.notify();
            }
            _ => {}
        });

        Self {
            graph,
            _subscription: subscription,
            build_task: None,
            model: None,
            snapshot: None,
            error: None,
            loading: false,
            limit: 250,
            selected: 0,
            previous: None,
            edge: None,
            context_tab: false,
            expanded_field: None,
            state_page: 0,
            neighborhood: false,
            edge_labels: false,
            fit_pending: false,
        }
    }

    fn open_model(&mut self, model: Model, cx: &mut Context<Self>) {
        self.model = Some(model);
        self.snapshot = None;
        self.error = None;
        self.loading = true;
        self.selected = 0;
        self.previous = None;
        self.edge = None;
        self.state_page = 0;
        self.expanded_field = None;
        self.graph.update(cx, |g, cx| g.clear(cx));
        let limit = self.limit;
        let work = cx
            .background_executor()
            .spawn(async move { model.build(limit) });
        self.build_task = Some(cx.spawn(async move |this, cx| {
            let result = work.await;
            let _ = this.update(cx, |this, cx| {
                this.loading = false;
                match result {
                    Ok(snapshot) => {
                        this.snapshot = Some(snapshot);
                        this.refresh_graph(cx);
                    }
                    Err(error) => this.error = Some(error),
                }
                cx.notify();
            });
        }));
        cx.notify();
    }

    fn dashboard(&mut self, cx: &mut Context<Self>) {
        self.build_task = None;
        self.loading = false;
        self.model = None;
        self.snapshot = None;
        self.error = None;
        self.graph.update(cx, |g, cx| g.clear(cx));
        cx.notify();
    }

    fn refresh_graph(&mut self, cx: &mut Context<Self>) {
        let Some(snapshot) = &self.snapshot else {
            return;
        };

        let mut visible = vec![!self.neighborhood; snapshot.states.len()];
        visible[self.selected] = true;
        if self.neighborhood {
            for edge in &snapshot.transitions {
                if edge.from == self.selected || edge.to == self.selected {
                    visible[edge.from] = true;
                    visible[edge.to] = true;
                }
            }
        }
        let nodes = snapshot
            .states
            .iter()
            .enumerate()
            .filter(|(i, _)| visible[*i])
            .map(|(i, s)| {
                Node::new(
                    i.to_string(),
                    format!("#{i} {}", s.label),
                    s.category.clone(),
                )
            })
            .collect();
        let edges = snapshot
            .transitions
            .iter()
            .enumerate()
            .filter(|(_, e)| visible[e.from] && visible[e.to])
            .map(|(i, e)| {
                Edge::new(e.from.to_string(), e.to.to_string(), e.step.clone())
                    .property("index", i.to_string())
            })
            .collect();
        let id = self.selected.to_string();
        self.graph.update(cx, |g, cx| {
            g.set_data(GraphData::new(nodes, edges), cx);
            g.select_node(Some(&id), false, cx);
            g.fit_view(cx);
        });
        self.fit_pending = true;
    }

    fn select(&mut self, index: usize, center: bool, cx: &mut Context<Self>) {
        if self
            .snapshot
            .as_ref()
            .is_none_or(|s| index >= s.states.len())
        {
            return;
        }

        let changed = self.selected != index;
        if changed {
            self.previous = Some(self.selected);
            self.selected = index;
            self.expanded_field = None;
        }
        self.edge = None;
        self.context_tab = false;
        self.state_page = index / PAGE_SIZE;
        if self.neighborhood && changed {
            self.refresh_graph(cx);
        } else if center {
            let id = index.to_string();
            self.graph
                .update(cx, |g, cx| g.select_node(Some(&id), true, cx));
        }
        cx.notify();
    }

    fn home(&self, cx: &mut Context<Self>) -> impl IntoElement {
        div()
            .id("dashboard")
            .flex_1()
            .overflow_y_scroll()
            .p_8()
            .flex()
            .flex_col()
            .gap_5()
            .child(div().text_3xl().child("Explore a system"))
            .child(
                div()
                    .text_color(rgb(MUTED))
                    .child("Each node is a state. Each directed edge is a valid transformation."),
            )
            .child(div().flex().flex_wrap().gap_4().children(
                Model::ALL.into_iter().enumerate().map(|(i, model)| {
                    button(format!("model-{i}"), "")
                        .accessibility_label(model.title())
                        .w(px(355.))
                        .min_h(px(150.))
                        .p_5()
                        .flex()
                        .flex_col()
                        .items_start()
                        .gap_3()
                        .child(div().text_lg().text_color(rgb(ACCENT)).child(model.title()))
                        .child(
                            div()
                                .text_sm()
                                .text_color(rgb(MUTED))
                                .child(model.description()),
                        )
                        .child(div().text_sm().child("Open graph →"))
                        .on_click(cx.listener(move |this, _, _, cx| this.open_model(model, cx)))
                }),
            ))
            .child(heading("Model sketches"))
            .child(div().text_sm().text_color(rgb(MUTED)).child(
                "These representations exist in the source but do not yet implement transitions.",
            ))
            .children(
                [
                    (
                        "Classical molecular dynamics",
                        "MolConformation · torsions / atoms and bonds",
                    ),
                    (
                        "Atoms & diatomic molecules",
                        "Atom, DiatomicMol · electron energy levels",
                    ),
                    (
                        "Spectroscopic states",
                        "SpectroState · term, vibration, rotation",
                    ),
                    (
                        "Density lattice",
                        "DensityLattice · spin density and nuclear sites",
                    ),
                    (
                        "Chemical bond graph",
                        "ChemGraph · atoms, bonds, and orbital occupations",
                    ),
                ]
                .into_iter()
                .map(|(title, detail)| {
                    div()
                        .flex()
                        .justify_between()
                        .p_3()
                        .border_b_1()
                        .border_color(rgb(LINE))
                        .child(
                            div()
                                .child(title)
                                .child(div().text_sm().text_color(rgb(MUTED)).child(detail)),
                        )
                        .child(
                            div()
                                .text_sm()
                                .text_color(rgb(MUTED))
                                .child("Not implemented"),
                        )
                }),
            )
    }

    fn sidebar(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let mut panel = div()
            .w(px(230.))
            .h_full()
            .flex_shrink_0()
            .bg(rgb(PANEL))
            .border_r_1()
            .border_color(rgb(LINE))
            .flex()
            .flex_col()
            .p_3()
            .gap_2()
            .child(
                button("home", "← All models")
                    .on_click(cx.listener(|this, _, _, cx| this.dashboard(cx))),
            )
            .child(heading("Systems"));
        for (i, model) in Model::ALL.into_iter().enumerate() {
            panel = panel.child(
                button(format!("switch-{i}"), model.title())
                    .bg(rgb(if self.model == Some(model) {
                        LINE
                    } else {
                        PANEL
                    }))
                    .on_click(cx.listener(move |this, _, _, cx| this.open_model(model, cx))),
            );
        }
        if let Some(snapshot) = &self.snapshot {
            let start = self.state_page * PAGE_SIZE;
            let end = (start + PAGE_SIZE).min(snapshot.states.len());
            panel = panel
                .child(heading(&format!(
                    "States · {}–{} of {}",
                    start + 1,
                    end,
                    snapshot.states.len()
                )))
                .child(
                    div()
                        .id("states")
                        .flex_1()
                        .min_h_0()
                        .overflow_y_scroll()
                        .flex()
                        .flex_col()
                        .gap_1()
                        .children(snapshot.states[start..end].iter().enumerate().map(
                            |(offset, state)| {
                                let i = start + offset;

                                button(format!("state-{i}"), format!("#{i} {}", state.label))
                                    .text_xs()
                                    .bg(rgb(if self.selected == i { LINE } else { PANEL }))
                                    .on_click(
                                        cx.listener(move |this, _, _, cx| this.select(i, true, cx)),
                                    )
                            },
                        )),
                )
                .child(
                    div()
                        .flex()
                        .gap_2()
                        .child(
                            button("prev-page", "Previous")
                                .disabled(start == 0)
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.state_page = this.state_page.saturating_sub(1);
                                    cx.notify();
                                })),
                        )
                        .child(
                            button("next-page", "Next")
                                .disabled(end == snapshot.states.len())
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.state_page += 1;
                                    cx.notify();
                                })),
                        ),
                );
        }

        panel
    }

    fn toolbar(&self, cx: &mut Context<Self>) -> impl IntoElement {
        div()
            .flex()
            .flex_wrap()
            .gap_2()
            .p_3()
            .border_b_1()
            .border_color(rgb(LINE))
            .child(button("fit", "Fit").on_click(
                cx.listener(|this, _, _, cx| this.graph.update(cx, |g, cx| g.fit_view(cx))),
            ))
            .child(
                button("zoom-out", "−")
                    .accessibility_label("Zoom out")
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.graph.update(cx, |g, cx| g.zoom_by(0.8, cx))
                    })),
            )
            .child(
                button("zoom-in", "+")
                    .accessibility_label("Zoom in")
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.graph.update(cx, |g, cx| g.zoom_by(1.25, cx))
                    })),
            )
            .child(button("relayout", "Arrange").on_click(
                cx.listener(|this, _, _, cx| this.graph.update(cx, |g, cx| g.relayout(cx))),
            ))
            .child(
                button(
                    "local",
                    if self.neighborhood {
                        "Neighborhood ✓"
                    } else {
                        "Neighborhood"
                    },
                )
                .on_click(cx.listener(|this, _, _, cx| {
                    this.neighborhood = !this.neighborhood;
                    this.refresh_graph(cx);
                    cx.notify();
                })),
            )
            .child(
                button(
                    "edge-labels",
                    if self.edge_labels {
                        "Steps ✓"
                    } else {
                        "Steps"
                    },
                )
                .on_click(cx.listener(|this, _, _, cx| {
                    this.edge_labels = !this.edge_labels;
                    this.graph
                        .update(cx, |g, cx| g.set_edge_labels(this.edge_labels, cx));
                    cx.notify();
                })),
            )
            .child(
                button("limit", format!("Limit: {}", self.limit))
                    .disabled(self.loading)
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.limit = match this.limit {
                            250 => 1000,
                            1000 => 5000,
                            _ => 250,
                        };
                        if this.limit == 5000 {
                            this.neighborhood = true;
                        }
                        if let Some(model) = this.model {
                            this.open_model(model, cx);
                        }
                    })),
            )
    }

    fn field(
        &self,
        section: &str,
        field: &Field,
        before: Option<&Field>,
        cx: &mut Context<Self>,
    ) -> gpui::Div {
        let key = format!("{section}-{}", field.name);
        let expanded = self.expanded_field.as_deref() == Some(&key);
        let changed = before.is_some_and(|old| old.value != field.value);
        let preview = compact(&field.value, 100);
        let toggle_key = key.clone();
        let copy = field.value.clone();
        let mut row = div()
            .flex()
            .flex_col()
            .gap_1()
            .py_2()
            .border_b_1()
            .border_color(rgb(LINE))
            .child(
                button(key.clone(), "")
                    .accessibility_label(field.name.clone())
                    .w_full()
                    .flex()
                    .flex_col()
                    .items_start()
                    .gap_1()
                    .bg(rgb(PANEL))
                    .child(
                        div()
                            .text_xs()
                            .text_color(rgb(if changed { CHANGED } else { MUTED }))
                            .child(format!(
                                "{} {}",
                                if expanded { "▾" } else { "▸" },
                                field.name
                            )),
                    )
                    .child(
                        div()
                            .text_sm()
                            .text_color(rgb(if changed { CHANGED } else { TEXT }))
                            .child(preview),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.expanded_field = if this.expanded_field.as_deref() == Some(&toggle_key)
                        {
                            None
                        } else {
                            Some(toggle_key.clone())
                        };
                        cx.notify();
                    })),
            );
        if changed {
            row = row.child(
                div()
                    .px_3()
                    .text_xs()
                    .text_color(rgb(CHANGED))
                    .child(format!("was {}", compact(&before.unwrap().value, 100))),
            );
        }
        if expanded {
            row = row
                .child(
                    div()
                        .id(SharedString::from(format!("value-{key}")))
                        .max_h(px(260.))
                        .overflow_scroll()
                        .p_3()
                        .bg(rgb(BG))
                        .child(
                            div()
                                .font_family("Consolas")
                                .text_xs()
                                .child(field.value.clone()),
                        ),
                )
                .child(
                    button(format!("copy-{key}"), "Copy value").on_click(move |_, _, cx| {
                        cx.write_to_clipboard(ClipboardItem::new_string(copy.clone()))
                    }),
                );
        }

        row
    }

    fn inspector(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let mut panel = div()
            .w(px(350.))
            .h_full()
            .flex_shrink_0()
            .bg(rgb(PANEL))
            .border_l_1()
            .border_color(rgb(LINE))
            .flex()
            .flex_col()
            .p_4()
            .gap_3()
            .child(
                div()
                    .flex()
                    .gap_2()
                    .child(
                        button("state-tab", "State")
                            .bg(rgb(if self.context_tab { PANEL } else { LINE }))
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.context_tab = false;
                                this.expanded_field = None;
                                cx.notify();
                            })),
                    )
                    .child(
                        button("context-tab", "Context")
                            .bg(rgb(if self.context_tab { LINE } else { PANEL }))
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.context_tab = true;
                                this.expanded_field = None;
                                cx.notify();
                            })),
                    ),
            );
        let Some(snapshot) = &self.snapshot else {
            return panel.child(
                div()
                    .text_color(rgb(MUTED))
                    .child("Select a system to inspect its fields."),
            );
        };

        let mut body = div()
            .id("inspector-body")
            .flex_1()
            .min_h_0()
            .overflow_y_scroll()
            .flex()
            .flex_col()
            .gap_2();
        if self.context_tab {
            body = body
                .child(heading("Static context"))
                .child(
                    div()
                        .text_xs()
                        .text_color(rgb(MUTED))
                        .child("Fixed for this graph. Click a field to expand or copy it."),
                )
                .children(
                    snapshot
                        .context
                        .iter()
                        .map(|f| self.field("context", f, None, cx)),
                );
        } else {
            let (index, compare) = self.edge.map_or((self.selected, self.previous), |i| {
                let edge = &snapshot.transitions[i];

                (edge.to, Some(edge.from))
            });
            let state = &snapshot.states[index];
            body = body.child(heading(&format!("State #{index}")));
            if let Some(i) = self.edge {
                let edge = &snapshot.transitions[i];
                let (from, to) = (edge.from, edge.to);
                body = body
                    .child(
                        div()
                            .text_sm()
                            .text_color(rgb(ACCENT))
                            .child(format!("#{from} → #{to}\n{}", edge.step)),
                    )
                    .child(
                        div()
                            .flex()
                            .gap_2()
                            .child(button("source", "Go to source").on_click(
                                cx.listener(move |this, _, _, cx| this.select(from, true, cx)),
                            ))
                            .child(button("target", "Go to target").on_click(
                                cx.listener(move |this, _, _, cx| this.select(to, true, cx)),
                            )),
                    );
            }
            if let Some(previous) = compare {
                body = body.child(
                    div()
                        .text_xs()
                        .text_color(rgb(CHANGED))
                        .child(format!("Changed fields compared with state #{previous}")),
                );
            }
            body = body
                .child(heading("Transforming fields"))
                .children(state.fields.iter().map(|f| {
                    let before = compare.and_then(|i| {
                        snapshot.states[i]
                            .fields
                            .iter()
                            .find(|old| old.name == f.name)
                    });

                    self.field("state", f, before, cx)
                }));
            if !state.observables.is_empty() {
                body =
                    body.child(heading("Computed values"))
                        .children(state.observables.iter().map(|f| {
                            let before = compare.and_then(|i| {
                                snapshot.states[i]
                                    .observables
                                    .iter()
                                    .find(|old| old.name == f.name)
                            });

                            self.field("derived", f, before, cx)
                        }));
            }
            let outgoing = snapshot
                .transitions
                .iter()
                .enumerate()
                .filter(|(_, e)| e.from == index)
                .collect::<Vec<_>>();
            body = body
                .child(heading(&format!("Outgoing steps · {}", outgoing.len())))
                .child(div().text_xs().text_color(rgb(MUTED)).child(
                    "Click a step to follow it. Click a graph edge to compare its endpoints.",
                ));
            if outgoing.is_empty() {
                body = body.child(div().text_sm().text_color(rgb(MUTED)).child(
                    if snapshot.truncated {
                        "No outgoing steps in this partial graph. Raise the limit to explore more."
                    } else {
                        "No valid outgoing steps."
                    },
                ));
            }
            for (i, edge) in outgoing {
                let to = edge.to;
                body = body.child(
                    button(format!("step-{i}"), format!("→ #{to}  {}", edge.step))
                        .text_xs()
                        .on_click(cx.listener(move |this, _, _, cx| {
                            this.select(index, false, cx);
                            this.select(to, true, cx);
                        })),
                );
            }
        }
        panel = panel.child(body);

        panel
    }
}

impl Render for Explorer {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        if std::mem::take(&mut self.fit_pending) {
            let graph = self.graph.clone();
            window.on_next_frame(move |_, cx| graph.update(cx, |g, cx| g.fit_view(cx)));
        }

        let mut root = div()
            .size_full()
            .flex()
            .flex_col()
            .bg(rgb(BG))
            .text_color(rgb(TEXT))
            .font_family("Segoe UI")
            .child(
                div()
                    .h(px(62.))
                    .flex_shrink_0()
                    .px_6()
                    .flex()
                    .items_center()
                    .justify_between()
                    .border_b_1()
                    .border_color(rgb(LINE))
                    .child(div().text_lg().child("STATE / GRAPH"))
                    .child(
                        div()
                            .text_sm()
                            .text_color(rgb(MUTED))
                            .child(self.model.map_or("Model dashboard", Model::title)),
                    ),
            );
        if self.model.is_none() {
            return root.child(self.home(cx));
        }

        let mut center = div()
            .flex_1()
            .min_w_0()
            .h_full()
            .flex()
            .flex_col()
            .child(self.toolbar(cx));
        if self.loading {
            center = center.child(
                div()
                    .flex_1()
                    .flex()
                    .flex_col()
                    .items_center()
                    .justify_center()
                    .gap_3()
                    .child("Generating reachable states…")
                    .child(
                        div()
                            .text_sm()
                            .text_color(rgb(MUTED))
                            .child("You can switch models while this runs."),
                    ),
            );
        } else if let Some(error) = &self.error {
            center = center.child(
                div()
                    .p_6()
                    .child(format!("Could not generate this graph: {error}"))
                    .child(
                        button("retry", "Retry").on_click(cx.listener(|this, _, _, cx| {
                            if let Some(model) = this.model {
                                this.open_model(model, cx);
                            }
                        })),
                    ),
            );
        } else {
            center = center.child(
                div()
                    .flex_1()
                    .min_h_0()
                    .overflow_hidden()
                    .child(self.graph.clone()),
            );
        }
        if let Some(snapshot) = &self.snapshot {
            let categories = self.graph.read(cx).categories();
            let status = if snapshot.truncated {
                "Partial graph · generation limit reached"
            } else {
                "Complete reachable graph"
            };
            center = center.child(div().p_3().text_xs().border_t_1().border_color(rgb(LINE)).flex().flex_col().gap_1()
                .child(div().flex().flex_wrap().gap_3().children(categories.into_iter().map(|category| {
                    div().text_color(category.color).child(format!("● {} · {}", category.name, category.count))
                })))
                .child(div().text_color(rgb(if snapshot.truncated { CHANGED } else { ACCENT }))
                    .child(format!("{} states · {} transitions · {status}", snapshot.states.len(), snapshot.transitions.len())))
                .child(div().text_color(rgb(MUTED)).child("Drag to pan · wheel / Ctrl+scroll to zoom · drag nodes to arrange · arrows to select · 0 to fit")));
        }
        root = root.child(
            div()
                .flex_1()
                .min_h_0()
                .flex()
                .child(self.sidebar(cx))
                .child(center)
                .child(self.inspector(cx)),
        );

        root
    }
}

fn button(id: impl Into<SharedString>, label: impl Into<SharedString>) -> Button {
    let label = label.into();

    Button::new(id.into())
        .accessibility_label(label.clone())
        .px_3()
        .py_2()
        .rounded_md()
        .bg(rgb(PANEL))
        .text_sm()
        .hover(|s| s.bg(rgb(LINE)))
        .focus_visible(|s| s.border_1().border_color(rgb(ACCENT)))
        .styles(|s| s.disabled(|s| s.opacity(0.4)))
        .when(!label.is_empty(), |button| button.child(label))
}

fn heading(text: &str) -> impl IntoElement {
    div()
        .mt_3()
        .mb_1()
        .text_sm()
        .text_color(rgb(MUTED))
        .child(text.to_owned())
}

fn compact(value: &str, max: usize) -> String {
    let flat = value.split_whitespace().collect::<Vec<_>>().join(" ");
    if flat.chars().count() <= max {
        return flat;
    }

    format!("{}…", flat.chars().take(max).collect::<String>())
}
