import React from "react";
import {Col, Container, Row, Badge, Placeholder, Accordion, Popover, OverlayTrigger} from "react-bootstrap";
import ApexChart from "../ApexChart/ApexChart";
import './SearchResultObject.css'
import DrugDescription from "../DrugDescription/DrugDescription";


const SearchResultObject = ( props ) => {
    const { termCountDict, totalCount } = props.searchResults.result
    const eventsOverTime = props.eventsOverTime.result.results
    const [retrievedTermArr, setRetrievedTermArr] = React.useState(null)

    let groupDescription = ''

    const processDrugGroups = (groups) => {
        return groups.map((group, index) => {
            let variant
            switch (group) {
                case 'approved':
                    variant = 'success'
                    groupDescription = 'Approved'
                    break
                case 'investigational':
                    variant = 'warning'
                    groupDescription = 'Investigational'
                    break
                case 'illicit':
                    variant = 'danger'
                    groupDescription = 'Illicit'
                    break
                case 'vet_approved':
                    variant = 'info'
                    groupDescription = 'Veterinary Approved'
                    break
                default:
                    variant = 'secondary'
            }
            group = group.substring(0, 3).toUpperCase()

            const popover = (
                <Popover id={`popover-${index}`}>
                    <Popover.Body>
                        {groupDescription}
                    </Popover.Body>
                </Popover>
            )

            return (
                <OverlayTrigger key={index} trigger={['hover', 'focus']} placement={'bottom'} overlay={popover}>
                    <Badge className={'mx-1 no-pointer-badge'} key={index} bg={variant}>{group}</Badge>
                </OverlayTrigger>
            )
        })

    }

    return (
        <div>
            <Container className="my-4">
                {retrievedTermArr ? (
                    <>
                        {retrievedTermArr && <h1>{retrievedTermArr.drug_name}</h1>}
                        {retrievedTermArr && <h3>{processDrugGroups(retrievedTermArr.groups)}</h3>}
                    </>
                ) : (
                    <>
                        <Placeholder as={'h1'} animation="glow">
                            <Placeholder xs={6} />
                        </Placeholder>
                        <Placeholder as={'h3'} animation="glow">
                            <Placeholder xs={3} />
                        </Placeholder>
                    </>
                )}
            </Container>
            <Container className={'main-drug-container'}>
                <Row>
                    <Col className={'alert alert-secondary mx-2'}>
                        <span>The provided information is indicative and shouldn't be used for inference.</span>
                    </Col>
                </Row>
                <Row>
                    <Col>
                        <DrugDescription drugName={props.searchTerm} onRetrieved={setRetrievedTermArr}/>
                    </Col>
                    <Col>
                        <Accordion>
                            <Accordion.Item eventKey={'0'}>
                                <Accordion.Header><Badge>ADEs Reported</Badge></Accordion.Header>
                                <Accordion.Body>
                                    {totalCount}
                                </Accordion.Body>
                            </Accordion.Item>
                            <Accordion.Item eventKey={'1'}>
                                <Accordion.Header><Badge>Indication</Badge></Accordion.Header>
                                <Accordion.Body>
                                    {retrievedTermArr ? (
                                        retrievedTermArr.indication.split(/(?<=\.)/)[0]
                                    ) : (
                                        <Placeholder animation="glow">
                                            <Placeholder xs={7} /> <Placeholder xs={4} /> <Placeholder xs={4} />{' '}
                                            <Placeholder xs={6} /> <Placeholder xs={8} />
                                        </Placeholder>
                                    )}
                                </Accordion.Body>
                            </Accordion.Item>
                            <Accordion.Item eventKey={'2'}>
                                <Accordion.Header><Badge>Half-life</Badge></Accordion.Header>
                                <Accordion.Body>
                                    {retrievedTermArr ? (
                                        retrievedTermArr.half_life
                                    ) : (
                                        <Placeholder animation="glow">
                                            <Placeholder xs={7} /> <Placeholder xs={4} /> <Placeholder xs={4} />{' '}
                                            <Placeholder xs={6} />
                                        </Placeholder>
                                    )}
                                </Accordion.Body>
                            </Accordion.Item>
                        </Accordion>
                    </Col>
                </Row>
                <Row>
                    <Col className={'mt-3'}>
                        {props.eventsOverTime && <ApexChart eventDict={eventsOverTime} type={'events_over_time'} />}
                        {props.searchResults && <ApexChart eventDict={termCountDict} totalCount={totalCount} type={'all_groups'} />}
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default SearchResultObject;